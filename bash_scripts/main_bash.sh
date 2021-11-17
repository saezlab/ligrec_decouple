#!/bin/bash
## The One bash script to rule them all

## Submit BRCA scripts
bash comp_brca_prep ER
bash comp_brca_prep TNBC
bash comp_brca_prep HER2

### Run Mouse Brain
bash large_bash.sh brain

## Submit CITE-Seq pipe
bash large_bash.sh analysis/citeseq/citeseq_prep.R citeseq

## Submit Comparisons
### Script to Run
compscript = "ligrec_decouple/analysis/comparison/liana_comparison_pipe.R"

### Run CRC
bash generic_bash.sh compscript crc data/input/comparison/crc_data/crc_korean_form.rds

### Run CBMC
bash generic_bash.sh compscript cbmc data/input/spatial/brain_cortex/brain_liana_results.RDS

### Run Pancreatic
bash generic_bash.sh compscript panc8 data/input/comparison/panc8.RDS

