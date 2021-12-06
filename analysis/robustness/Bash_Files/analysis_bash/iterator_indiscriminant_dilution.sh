#!/bin/bash
#SBATCH	-p single
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --mem=60000
#SBATCH --job-name="modify_iterator"
#SBATCH --output=modify_iterator.out
#SBATCH --mail-user=daniel.dimitrov@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/bq_ddimitrov/Repos/ligrec_decouple


# Add arguments to pass to Run_Iterator. 1 = modify_baseline, 2 = job_id
RScipt analysis/robustness/Code/Analysis_Scripts/Resources_Run_Iterator.R TRUE $SLURM_JOBID
