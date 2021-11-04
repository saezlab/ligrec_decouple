# Citeseq-LR Correlations and Specificity ROC

## Reproduce CiteSeq Analysis
```
Run citeseq_convert.Rmd - convert murine_datasets.h5a to seurat objects

# Run to perform basic analysis on the 10x datasets, LR-ADT_Receptor correlations,
# and (z-normalized) Receptor abundance specificity ROC
```

Generate Seurat objects and LAINA results
```
source("analysis/citeseq/citeseq_prep.R")
```

Generate Plots
```
source("analysis/citeseq/citeseq_analyse.R")
```


## Datasets

### 5k PBMCs from 10x
##### 5k PBMCs from a Healthy Donor (NextGem)
https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.1.0/5k_pbmc_protein_v3_nextgem
##### 5k PBMCs from a Healthy Donor
https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.1.0/5k_pbmc_protein_v3
*31 surface antibodies each

### 10k PBMCs and MALT
#####
https://www.10xgenomics.com/resources/datasets/10-k-pbm-cs-from-a-healthy-donor-gene-expression-and-cell-surface-protein-3-standard-3-0-0
##### 10k Cells from a MALT Tumor
https://www.10xgenomics.com/resources/datasets/10-k-cells-from-a-malt-tumor-gene-expression-and-cell-surface-protein-3-standard-3-0-0
*18 surface antibodies each


### scRNAseq and 13-antibody sequencing of CBMCs, acquired from [`SeuratData`v0.2.1]:
```{r}
SeuratData::InstallData("cbmc")
```

### Murine Spleen-Lymph CiteSeq data
2 datasets: one with 111 and another with 208 ADTs - `SLN111` and `SLN208`, respectively
Processed counts were obtained from:
https://www.nature.com/articles/s41592-020-01050-x#data-availability
