#!/bin/bash
## Bash to run all script for all comparison summaries
summscript="analysis/comparison/comparison_summarise.R"

## Define Variables (setting = "specs_frac", "house_n", "specs_n", 'comp_n', 'comp_frac')
setting=$1

### Run CRC
bash bash_scripts/generic_bash.sh $summscript "crc_"$setting "data/output/comparison_out/crc_liana_res.RDS" $setting

### Run CBMC
bash bash_scripts/generic_bash.sh $summscript "cbmc_"$setting "data/output/comparison_out/cbmc_liana_res.RDS" $setting

### Run Pancreatic
bash bash_scripts/generic_bash.sh $summscript "panc8_"$setting "data/output/comparison_out/panc8_liana_res.RDS" $setting

### Run ER (BRCA)
bash bash_scripts/generic_bash.sh $summscript "er_"$setting "data/output/comparison_out/BRCA_ER_liana_res.RDS" $setting

### Run HER2 (BRCA)
bash bash_scripts/generic_bash.sh $summscript "her2_"$setting "data/output/comparison_out/BRCA_HER2_liana_res.RDS" $setting

### Run TNBC (BRCA)
bash bash_scripts/generic_bash.sh $summscript "tnbc_"$setting "data/output/comparison_out/BRCA_TNBC_liana_res.RDS" $setting
