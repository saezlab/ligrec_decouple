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

# $SLURM_SUBMIT_DIR project directory, $1 cancer subtype  - subtypes=("ER" "TNBC" "HER2")
 echo "Submitted: Rscript analysis/comparison/comp_brca_prep.R $SLURM_SUBMIT_DIR $1"
 Rscript analysis/comparison/comp_brca_prep.R $SLURM_SUBMIT_DIR $1
