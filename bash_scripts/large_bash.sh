#!/bin/bash
sbatch <<EOT
#!/bin/bash

#SBATCH -p single
#SBATCH -N 1
#SBATCH --time=71:00:00
#SBATCH --mem=200000
#SBATCH --job-name=$2"_brca_run"
#SBATCH --output="out/"$2".out"
#SBATCH --error="out/"$2".err"
#SBATCH --mail-user=daniel.dimitrov@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/bq_ddimitrov/Repos/ligrec_decouple

echo "Submitted: $1 (Job ID: $SLURM_JOBID)"
Rscript $1

hostname
exit 0

EOT
