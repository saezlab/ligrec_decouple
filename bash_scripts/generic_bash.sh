#!/bin/bash
sbatch <<EOT
#!/bin/bash

#SBATCH -p single
#SBATCH -N 1
#SBATCH --time=47:00:00
#SBATCH --mem=40000
#SBATCH --job-name=$1"_brca_run"
#SBATCH --output="out/"$1".out"
#SBATCH --error="out/"$1".err"
#SBATCH --mail-user=daniel.dimitrov@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/bq_ddimitrov/Repos/ligrec_decouple

echo "Submitted: $1 (Job ID: $SLURM_JOBID)"
Rscript $1

hostname
exit 0

EOT
