#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH --time=23:00:00
#SBATCH --mem=40000
#SBATCH --output=out/R-%x.%j.out
#SBATCH --mail-user=daniel.dimitrov@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/bq_ddimitrov/Repos/ligrec_decouple
#SBATCH --cpus-per-task 12

echo "Submitted: $1 (Job ID: $SLURM_JOBID)"
Rscript $1
