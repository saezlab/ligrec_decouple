#!/bin/bash
sbatch <<EOT
#!/bin/bash

#SBATCH -p single
#SBATCH -N 1
#SBATCH --time=71:00:00
#SBATCH --mem=100000
#SBATCH --job-name=$1"brca_run"
#SBATCH --output="out/"$1"brca_run.out"
#SBATCH --e="out/"$1"brca_run.err"
#SBATCH --mail-user=daniel.dimitrov@uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/bq_ddimitrov/Repos/ligrec_decouple

# project directory, $1 cancer subtype  - subtypes=("ER" "TNBC" "HER2")
echo "Submitted: Rscript analysis/comparison/comp_brca_prep.R $1"
Rscript analysis/comparison/comp_brca_prep.R $1

hostname
exit 0

EOT
