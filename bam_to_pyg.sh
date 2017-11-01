#!/bin/bash

# minimum overlap of 51% may cause slight differences with GOTHiC
bedtools pairtobed \
    -f 0.51 \
    -bedpe \
    -abam $1 \
    -b <( tail -n +3 $3 | cut -f1-3 ) \
    -type both \
    | cut -f11-13 \
    | paste - - \
    | gzip \
    > $2
