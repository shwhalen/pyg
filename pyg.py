#!/usr/bin/env python

import itertools
import json
import numpy as np
import pandas as pd
import sklearn.externals.joblib as joblib
import sys

from tqdm import tqdm

from rpy2.robjects import r, numpy2ri, pandas2ri
numpy2ri.activate()
pandas2ri.activate()

# scipy's stats.binom.sf disagrees with r for low counts
def get_pvalues(parameters):
    return [
        (int(count), probability, r['binom.test'](count - 1, total_read_counts, probability, alternative = 'greater').rx2('p.value')[0])
        for count, probability in parameters
    ]


config = json.load(open(sys.argv[1], 'rt'))
input_fn = config['input_fn']
output_fn = config['output_fn']
cistrans = config['cistrans']
bin1_size = config['bin1_size']
bin2_size = config['bin2_size']
self_ligation_distance = config['self_ligation_distance']
n_jobs = config['n_jobs']
debug = True

re1_columns = ['re1_chrom', 're1_start', 're1_end']
re2_columns = ['re2_chrom', 're2_start', 're2_end']
bin1_columns = ['bin1_chrom', 'bin1_locus']
bin2_columns = ['bin2_chrom', 'bin2_locus']

print('parsing reads pre-mapped to restriction fragments')
reads_df = pd.read_csv(
    input_fn,
    sep = '\t',
    header = None,
    names = re1_columns + re2_columns
)

if cistrans == 'cis':
    print('retaining only cis interactions')
    reads_df = reads_df.query('re1_chrom == re2_chrom')
elif cistrans == 'trans':
    print('retaining only trans interactions')
    reads_df = reads_df.query('re1_chrom != re2_chrom')

print('binning fragment reads')
reads_df['bin1_locus'] = ((reads_df['re1_start'] // bin1_size) * bin1_size).astype(int)
reads_df['bin2_locus'] = ((reads_df['re2_start'] // bin2_size) * bin2_size).astype(int)
reads_df.rename(columns = {'re1_chrom': 'bin1_chrom', 're2_chrom': 'bin2_chrom'}, inplace = True)

# required to match GOTHiC
print('filtering self-ligations')
reads_df = reads_df.query('(bin1_chrom != bin2_chrom) or (abs(bin2_locus - bin1_locus) > @self_ligation_distance)')

print('ensuring consistent bin ordering')
unordered_mask = reads_df.eval('(bin1_chrom > bin2_chrom) or (bin1_chrom == bin2_chrom and bin2_locus < bin1_locus)')
if unordered_mask.sum() != 0:
    unordered_counts_df = reads_df[unordered_mask]
    unordered_counts_df.columns = bin2_columns + bin1_columns
    ordered_counts_df = reads_df[~unordered_mask]
    reads_df = (
        pd.concat([ordered_counts_df, unordered_counts_df], ignore_index = True)
        [bin1_columns + bin2_columns + ['observed_count']]
    )
assert reads_df.eval('(bin1_chrom < bin2_chrom) or (bin1_chrom == bin2_chrom and bin2_locus >= bin1_locus)').all()

print('aggregate reads into counts')
counts_df = (
    reads_df
    .groupby(bin1_columns + bin2_columns)
    .size()
    .rename('observed_count')
    .reset_index()
)
del reads_df

print('creating bin names to simplify grouping/filtering')
counts_df['bin1_name'] = counts_df['bin1_chrom'].astype(str) + ':' + counts_df['bin1_locus'].astype(str)
counts_df['bin2_name'] = counts_df['bin2_chrom'].astype(str) + ':' + counts_df['bin2_locus'].astype(str)

print('removing same-bin interactions')
counts_df = counts_df.query('bin1_name != bin2_name')

print('computing absolute coverage')
bin1_coverage = counts_df.groupby('bin1_name')['observed_count'].sum()
bin2_coverage = counts_df.groupby('bin2_name')['observed_count'].sum()
absolute_coverage = (
    pd.concat([bin1_coverage, bin2_coverage], axis = 1)
    .fillna(0)
    .sum(axis = 1)
    .astype(int)
    .rename('absolute_coverage')
)

print('computing relative coverage')
relative_coverage_df = (
    absolute_coverage
    .div(absolute_coverage.sum())
    .rename('relative_coverage')
    .to_frame()
)

print('merging relative coverage')
counts_df = (
    pd.merge(counts_df, relative_coverage_df, left_on = 'bin1_name', right_index = True)
    .rename(columns = {'relative_coverage': 'bin1_relative_coverage'})
)
counts_df = (
    pd.merge(counts_df, relative_coverage_df, left_on = 'bin2_name', right_index = True)
    .rename(columns = {'relative_coverage': 'bin2_relative_coverage'})
)

print('computing correction factor')
pairs_per_chrom = {}
for chrom in counts_df['bin1_chrom'].unique():
    unique_bin1_loci = set(counts_df.query('bin1_chrom == @chrom')['bin1_name'])
    unique_bin2_loci = set(counts_df.query('bin2_chrom == @chrom')['bin2_name'])
    unique_loci = unique_bin1_loci | unique_bin2_loci
    pairs_per_chrom[chrom] = len(unique_loci)
pairs_per_chrom = pd.Series(pairs_per_chrom)

total_read_counts = counts_df['observed_count'].sum()                                      # numberOfReadPairs
total_unique_bins = len(set(counts_df['bin1_name']) | set(counts_df['bin2_name']))         # len(all_bins)
total_interactions = total_unique_bins**2                                                  # numberOfAllInteractions
total_nonredundant_pairs = (total_unique_bins * (total_unique_bins - 1)) // 2              # upperhalfBinNumber
total_nonredundant_cis_pairs = ((pairs_per_chrom**2).sum() - total_unique_bins) // 2       # cisBinNumber
total_nonredundant_trans_pairs = total_nonredundant_pairs - total_nonredundant_cis_pairs   # transBinNumber

if cistrans == 'all':
    diagonal_probability = (relative_coverage_df.squeeze()**2).sum()                       # diagonalProb
    correction_factor = 1 / (1 - diagonal_probability)                                     # probabilityCorrection
    padjust_n = total_nonredundant_pairs + total_unique_bins
    if debug:
        print(diagonal_probability)
elif cistrans == 'cis':
    correction_factor = total_nonredundant_pairs / total_nonredundant_cis_pairs
    padjust_n = total_nonredundant_cis_pairs + total_unique_bins
elif cistrans == 'trans':
    correction_factor = total_nonredundant_pairs / total_nonredundant_trans_pairs
    padjust_n = total_nonredundant_trans_pairs

if debug:
    print(total_read_counts, total_unique_bins, total_interactions, total_nonredundant_pairs, total_nonredundant_cis_pairs, total_nonredundant_trans_pairs, correction_factor)

print('computing probability of loci randomly interacting based on relative coverage')
counts_df.eval('random_interaction_probability = bin1_relative_coverage * bin2_relative_coverage * 2 * @correction_factor', inplace = True)

print('computing expected counts')
counts_df['expected_count'] = counts_df['random_interaction_probability'] * total_read_counts

# pre-compute value for all observed (count, probability) pairs so identical parameters are not re-computed (~20x speedup for large datasets)
# could exchange accuracy for time by trimming probabilities to fixed # digits so fewer unique values
print('computing significance of interactions using binomial test')
parameters = counts_df[['observed_count', 'random_interaction_probability']].drop_duplicates().values
parameter_chunks = np.array_split(parameters, min(n_jobs, len(parameters)))
results = joblib.Parallel(n_jobs)(
    joblib.delayed(get_pvalues)(_)
    for _ in tqdm(parameter_chunks, 'parameter chunk')
)

print('merging pvalues')
pvalues_df = pd.DataFrame(
    list(itertools.chain.from_iterable(results)),
    columns = ['observed_count', 'random_interaction_probability', 'pvalue']
)
counts_df = pd.merge(counts_df, pvalues_df, on = ['observed_count', 'random_interaction_probability'])

print('correcting for multiple testing (n = {})'.format(padjust_n))
counts_df['qvalue'] = r['p.adjust'](counts_df['pvalue'], method = 'BH', n = padjust_n)

print('saving')
(
    counts_df
    .drop(['bin1_name', 'bin2_name'], axis = 1)
    .sort_values(['bin1_chrom', 'bin1_locus', 'bin2_chrom', 'bin2_locus'])
    .reset_index(drop = True)
    .to_feather(output_fn)
)
