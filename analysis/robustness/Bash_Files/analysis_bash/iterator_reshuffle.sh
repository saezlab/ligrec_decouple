#!/bin/bash
#SBATCH	-p single
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --mem=20000
#SBATCH --job-name="reshuffle_iterator"
#SBATCH --output=reshuffle_iterator.out
#SBATCH --mail-user=[insert email]
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/[insert filepath]/repos/ligrec_decouple


# Add arguments to pass to Run_Iterator. 1 = reshuffle_or_subset, 2 = job_id
/net/data.isilon/ag-saez/[filepath to liana_env and Rscript] analysis/robustness/Code/Analysis_Scripts/Clusters_Run_Iterator.R "reshuffle" $SLURM_JOBID  
