#!/bin/bash
#SBATCH -p single
#SBATCH -N 1
#SBATCH --time=71:00:00
#SBATCH --mem=200000
#SBATCH --output=out/R-%x.%j.out
#SBATCH --mail-user=daniel.dimitrov@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/bq_ddimitrov/Repos/ligrec_decouple

# Job name is passed with --job-name=$A.$b.run

echo "Submitted: $1 (Job ID: $SLURM_JOBID)"
Rscript $1
