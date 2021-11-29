#!/bin/bash
#SBATCH	-p single
#SBATCH -N 1
#SBATCH --time=1:00:00
#SBATCH --mem=8000
#SBATCH --job-name="seurat_tut"
#SBATCH --output=seurat_tut.out
#SBATCH --mail-user=p.burmedi@stud.uni-heidelberg.de
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --chdir /net/data.isilon/ag-saez/bq_pburmedi/repos/ligrec_robustness


/net/data.isilon/ag-saez/bq_pburmedi/SOFTWARE/miniconda3/envs/liana_env/bin/Rscript Code/Seurat_Tut_Script.R
