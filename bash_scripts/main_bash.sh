#!/bin/bash
## The One bash script to rule them all

## Submit BRCA scripts
bash comp_brca_prep ER
bash comp_brca_prep TNBC
bash comp_brca_prep HER2

## Submit Comparisons
### Run CRC
bash large_bash.sh

### Run CBMC
bash generic_bash.sh

### Run Pancreatic
bash generic_bash.sh

## Submit CITE-Seq pipe
bash large_bash.sh analysis/citeseq/citeseq_prep.R
