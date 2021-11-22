#!/bin/bash
## Bash to run all script for all comparison summaries
$summscript="analysis/comparison/comparison_analyse.R"

### Run CRC
# bash bash_scripts/large_bash.sh $summscript crc_out

### Run CBMC
# bash bash_scripts/large_bash.sh $summscript cbmc_out

### Run Pancreatic
bash bash_scripts/generic.sh $summscript panc8_sum "data/output/comparison_out/panc8_liana_res.RDS"

### Run ER (BRCA)
bash bash_scripts/large_bash.sh $summscript brca_er_sum "data/output/comparison_out/BRCA_ER_liana_res.RDS"

### Run HER2 (BRCA)
bash bash_scripts/large_bash.sh $summscript brca_her2_sum "data/output/comparison_out/BRCA_HER2_liana_res.RDS"

### Run TNBC (BRCA)
bash bash_scripts/large_bash.sh $summscript tnbc_tnbc_sum "data/output/comparison_out/BRCA_TNBC_liana_res.RDS"

