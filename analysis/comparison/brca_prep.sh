#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH --time=71:00:00
#SBATCH --mem=200000
#SBATCH --job-name="brca_run"
#SBATCH --output=out/brca_run.out
#SBATCH --mail-user=daniel.dimitrov@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/bq_ddimitrov/Repos/ligrec_decouple

## declare an array variable
declare -a subtypes=("ER" "TNBC" "HER2")

# get length of an array
arraylength=${#subtypes[@]}

# Loop over brca_types and first comm arg is path_to_project
for (( i=0; i<${arraylength}; i++ ));
do
 echo "Submitted: Rscript analysis/comparison/comp_brca_prep.R $1 ${subtypes[$i]}"
 Rscript analysis/comparison/comp_brca_prep.R $1 ${subtypes[$i]}
done
