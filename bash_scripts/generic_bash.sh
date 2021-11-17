#!/bin/bash
sbatch <<EOT
#!/bin/bash

#SBATCH -p single
#SBATCH -N 1
#SBATCH --time=23:00:00
#SBATCH --mem=40000
#SBATCH --job-name=$2
#SBATCH --output="out/"$2".out"
#SBATCH --error="out/"$2".err"
#SBATCH --mail-user=daniel.dimitrov@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/bq_ddimitrov/Repos/ligrec_decouple

## 1 Rscript 2 jobname 3 seurat_path
echo "Submitted: $1 (Job ID: $SLURM_JOBID; Job Name: $2)"
Rscript $1 $2 $3

hostname
exit 0

EOT
