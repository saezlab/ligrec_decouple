#!/bin/bash
#SBATCH	-p single
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --mem=60000
#SBATCH --job-name="baseline_iterator"
#SBATCH --output=baseline_iterator.out
#SBATCH --mail-user=[insert email]
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/[insert filepath]/repos/ligrec_decouple


# Add arguments to pass to Run_Iterator. 1 = modify_baseline, 2 = job_id
/net/data.isilon/ag-saez/[filepath to liana_env and Rscript] analysis/robustness/Code/Analysis_Scripts/Resources_Run_Iterator.R FALSE $SLURM_JOBID
