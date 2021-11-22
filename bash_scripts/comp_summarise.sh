#!/bin/bash
## Bash to run all script for all comparison summaries
summscript="analysis/comparison/comparison_summarise.R"

## Define Variables (setting = "specs_frac", "house_frac", "specs_n")
setting="specs_frac"

### Run CRC
## bash bash_scripts/large_bash.sh $summscript crc_out $setting

### Run CBMC
## bash bash_scripts/large_bash.sh $summscript cbmc_out $setting

### Run Pancreatic
bash bash_scripts/large_bash.sh $summscript "panc8_"$setting "data/output/comparison_out/panc8_liana_res.RDS" $setting

### Run ER (BRCA)
bash bash_scripts/large_bash.sh $summscript "er_"$setting "data/output/comparison_out/BRCA_ER_liana_res.RDS" $setting

### Run HER2 (BRCA)
bash bash_scripts/large_bash.sh $summscript "her2"$setting "data/output/comparison_out/BRCA_HER2_liana_res.RDS" $setting

### Run TNBC (BRCA)
bash bash_scripts/large_bash.sh $summscript "tnbc_"$setting "data/output/comparison_out/BRCA_TNBC_liana_res.RDS" $setting
