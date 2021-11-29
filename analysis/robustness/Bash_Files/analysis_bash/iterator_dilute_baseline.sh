#!/bin/bash
#SBATCH	-p single
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --mem=20000
#SBATCH --job-name="baseline_iterator"
#SBATCH --output=baseline_iterator.out
#SBATCH --mail-user=p.burmedi@stud.uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/bq_pburmedi/repos/ligrec_robustness


# Add arguments to pass to Run_Iterator. 1 = modify_baseline, 2 = feature_type, 3 = top_n, 4 = job_id, 5 = preserve_topology
/net/data.isilon/ag-saez/bq_pburmedi/SOFTWARE/miniconda3/envs/liana_env/bin/Rscript Code/Analysis_Scripts/Resources_Run_Iterator.R FALSE "variable" 500 $SLURM_JOBID FALSE 
