# LR correlation with ADTs

## Analysis
```
Run citeseq_convert.Rmd
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

*18 antibodies each


### scRNAseq and 13-antibody sequencing of CBMCs, acquired via:
```{r}
SeuratData::InstallData("cbmc")
```

### Public (murine) data
2 datasets: 1 with 111 antibodies and another with 208 ADTs (SLN111 and SLN208)
Available at:
https://www.nature.com/articles/s41592-020-01050-x#data-availability