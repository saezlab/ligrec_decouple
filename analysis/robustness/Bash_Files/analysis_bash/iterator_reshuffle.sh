#!/bin/bash
#SBATCH	-p single
#SBATCH -N 1
#SBATCH --time=72:00:00
#SBATCH --mem=60000
#SBATCH --job-name="reshuffle_iterator"
#SBATCH --output="out/reshuffle_iterator.out"
#SBATCH --mail-user=daniel.dimitrov@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/bq_ddimitrov/Repos/ligrec_decouple

# Add arguments to pass to Run_Iterator. 1 = reshuffle_or_subset, 2 = job_id
Rscript analysis/robustness/Code/Analysis_Scripts/Clusters_Run_Iterator.R "reshuffle" $SLURM_JOBID
