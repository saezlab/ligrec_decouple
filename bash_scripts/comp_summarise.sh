#!/bin/bash
## Bash to run all script for all comparison summaries
summscript="analysis/comparison/comparison_summarise.R"

## Define Variables (1 is top_frac/top_n; 2 is 0.05 or 1000)
top_fun=$1
top_x=$2

### Run CRC
## bash bash_scripts/large_bash.sh $summscript crc_out $top_fun $top_x

### Run CBMC
## bash bash_scripts/large_bash.sh $summscript cbmc_out $top_fun $top_x

### Run Pancreatic
bash bash_scripts/large_bash.sh $summscript panc8_sum "data/output/comparison_out/panc8_liana_res.RDS" $top_fun $top_x

### Run ER (BRCA)
bash bash_scripts/large_bash.sh $summscript brca_er_sum "data/output/comparison_out/BRCA_ER_liana_res.RDS" $top_fun $top_x

### Run HER2 (BRCA)
bash bash_scripts/large_bash.sh $summscript brca_her2_sum "data/output/comparison_out/BRCA_HER2_liana_res.RDS" $top_fun $top_x

### Run TNBC (BRCA)
bash bash_scripts/large_bash.sh $summscript tnbc_tnbc_sum "data/output/comparison_out/BRCA_TNBC_liana_res.RDS" $top_fun $top_x

