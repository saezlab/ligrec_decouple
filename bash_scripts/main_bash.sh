#!/bin/bash
## The One bash script to rule them all

## Submit BRCA scripts
bash bash_scripts/comp_brca_prep.sh ER
bash bash_scripts/comp_brca_prep.sh TNBC
bash bash_scripts/comp_brca_prep.sh HER2

### Run Mouse Brain
bash bash_scripts/large_bash.sh mouse_brain_prep_liana.R brain

## Submit CITE-Seq pipe
bash bash_scripts/large_bash.sh analysis/citeseq/citeseq_prep.R citeseq

## Submit Comparisons
### Script to Run Comparisons
compscript = "ligrec_decouple/analysis/comparison/comparisons_rscript.R"

### Run CRC
bash bash_scripts/generic_bash.sh compscript crc data/input/comparison/crc_data/crc_korean_form.rds

### Run CBMC
bash bash_scripts/generic_bash.sh compscript cbmc data/input/citeseq/cmbcs/cbmc_seurat.RDS

### Run Pancreatic
bash bash_scripts/generic_bash.sh compscript panc8 data/input/comparison/panc8.RDS

