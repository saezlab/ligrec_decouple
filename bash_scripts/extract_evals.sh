#!/bin/bash

extractscript="analysis/comparison/extract_aggregates.R"

## BRCA
### Extract aggregates Intersect
bash bash_scripts/generic_bash.sh $extractscript ER "data/output/comparison_out/BRCA_ER_liana_res.RDS" independent
bash bash_scripts/generic_bash.sh $extractscript HER2 "data/output/comparison_out/BRCA_HER2_liana_res.RDS" independent
bash bash_scripts/generic_bash.sh $extractscript TNBC "data/output/comparison_out/BRCA_TNBC_liana_res.RDS" independent

### Extract aggregates Intersect
bash bash_scripts/generic_bash.sh $extractscript ER "data/output/comparison_out/BRCA_ER_liana_res.RDS" intersect
bash bash_scripts/generic_bash.sh $extractscript HER2 "data/output/comparison_out/BRCA_HER2_liana_res.RDS" intersect
bash bash_scripts/generic_bash.sh $extractscript TNBC "data/output/comparison_out/BRCA_TNBC_liana_res.RDS" intersect

### Extract aggregates MAX
bash bash_scripts/generic_bash.sh $extractscript ER "data/output/comparison_out/BRCA_ER_liana_res.RDS" max
bash bash_scripts/generic_bash.sh $extractscript HER2 "data/output/comparison_out/BRCA_HER2_liana_res.RDS" max
bash bash_scripts/generic_bash.sh $extractscript TNBC "data/output/comparison_out/BRCA_TNBC_liana_res.RDS" max

## Brain
bash bash_scripts/generic_bash.sh $extractscript brain "data/input/spatial/brain_cortex/brain_liana_results.RDS" independent
bash bash_scripts/generic_bash.sh $extractscript brain "data/input/spatial/brain_cortex/brain_liana_results.RDS" intersect
bash bash_scripts/generic_bash.sh $extractscript brain "data/input/spatial/brain_cortex/brain_liana_results.RDS" max

### Run Cytosig Eval
# bash bash_scripts/large_bash.sh analysis/cytosig/cytosig_analyse.R cytosig

### Run Spatial Analyse
bash bash_scripts/large_bash.sh analysis/spatial/spatial_analyse.R space
