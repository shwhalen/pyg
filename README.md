# PyG

Python GOTHiC (PyG) is a Python implementation of the [GOTHiC](https://www.bioconductor.org/packages/release/bioc/html/GOTHiC.html) algorithm (v1.14) for assessing the statistical significance of chromatin interactions in Hi-C data.  We selected GOTHiC for its simplicity and use in many high profile papers.

This code was developed for the paper "The Epstein Barr virus episome maneuvers between nuclear 2 chromatin compartments during reactivation", and is published here for reproducibility.  Our analyses required genomic bins of different size in order to achieve better statistical power at informative resolutions, which the R implementation does not support, and thus necessitated this custom implementation.  In addition, we use `bedtools` to map reads to restriction fragments and compute the binomial test in parallel for non-redundant parameters, resulting in several orders of magnitude speedup and smaller memory footprint than the R implementation.

Unit tests show 99-100% correlation (Kendall's tau) between corresponding columns in the output files of GOTHiC and PyG.

## Dependencies

PyG primarily relies on `pandas`, `numpy`, and `bedtools`.  The `joblib` package (included in `scikit-learn`) is used for parallelization.  Differences between `scipy` and R's binomial test were observed; to ensure agreement with GOTHiC, we use the `rpy2` package for directly calling R's `binom.test` and `p.adjust` functions.

The dependencies for PyG, along with the versions used for testing, are:

- Python (3.6.3)
- R (3.4.2)
- pandas (0.20.3)
- numexpr (2.6.2)
- numpy (1.13.3)
- feather (0.4)
- scikit-learn (0.19.1)
- rpy2 (2.9)
- bedtools (2.26)

We recommend using [Miniconda](http://conda.pydata.org/miniconda.html) for R and Python packages.  [`bedtools`](https://github.com/arq5x/bedtools2/releases) has source and binary packages available for multiple platforms.

## Running

As with GOTHiC, we expect as input a BAM file of quality-controlled reads mapped to a genome produced by [HiCUP](https://www.bioinformatics.babraham.ac.uk/projects/hicup/).  Version 0.5.9 was used for the paper, paired with [bowtie2](https://github.com/BenLangmead/bowtie2/releases) 2.3.3.1 for alignment.  PyG maps HiCUP-processed reads to restriction fragments using the `bam_to_pyg.sh` script:

```
bam_to_pyg.sh bam_filename compressed_bed_output_filename digest_filename
```

where `bam_filename` is the BAM file produced by [HiCUP](https://www.bioinformatics.babraham.ac.uk/projects/hicup/) containing quality-controlled reads mapped to the genome, `compressed_output_filename` is the desired location of the output file (a GZIP'd bed file), and `digest_filename` is the name of the genome digest file produced by HiCUP's `hicup_digester` script.

PyG is run by passing a JSON configuration file:

```
pyg.py config.json
```

This file specifies all relevant parameters including input and output files, interaction types ("cis", "trans", or "all"), the size of interacting bins in # base pairs, a distance below which pairs will be considered self-ligated and thus filtered (10,000 base pairs is the default used by GOTHiC), and the number of CPU cores to use for the binomial test:

```
{
    "input_fn": "pyg_in/test_dataset1_2.pyg.gz",
    "output_fn": "pyg_out/test_dataset1_2.pyg.feather",
    "cistrans": "all",
    "bin1_size": 1e6,
    "bin2_size": 1e6,
    "self_ligation_distance": 10000,
    "n_jobs": 4
}
```

The output of PyG is a file in [feather](https://github.com/wesm/feather) format, which is a high performance container for columnar data with lightweight R and Python bindings.  This allows quick analysis of PyG output in either language.

## Unit Tests

The following are the correlation (Kendall's tau) and maximum difference between corresponding columns in GOTHiC and PyG output for the test Hi-C dataset included with HiCUP (mapped to hg19):

```
('relCoverage1', 'bin1_relative_coverage')                  
0.999302673593                
2.73188006363e-05             

('relCoverage2', 'bin2_relative_coverage')                  
0.998870240389                
2.73188006363e-05             

('probability', 'random_interaction_probability')           
0.997520954913                
4.1809365828e-08              

('readCount', 'observed_count')                             
1.0                           
0                             

('expected', 'expected_count')                              
0.997524800417                
0.000765229470417             

('pvalue', 'pvalue_')         
0.999980983819                
0.000755098835315             

('qvalue', 'qvalue_')         
1.0                           
0.00250407165144 
```

The following is the same for a larger dataset (GEO id SRX116341, mapped to mm10):

```
('relCoverage1', 'bin1_relative_coverage')                  
0.99943598081                 
3.13031425806e-07             

('relCoverage2', 'bin2_relative_coverage')                  
0.99959103775                 
3.13031425806e-07             

('probability', 'random_interaction_probability')           
0.999605553696                
5.73675839021e-10             

('readCount', 'observed_count')                             
0.999068702326                
30                            

('expected', 'expected_count')                              
0.999605553726                
0.0901615365302               

('pvalue', 'pvalue_')         
0.997679619907                
0.945555359362                

('qvalue', 'qvalue_')         
0.999663870462                
0.946159863912
```

Note the almost perfect correlation between columns, with a few differences due to using `bedtools` versus R for mapping reads to restriction fragments.