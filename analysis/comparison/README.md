## Comparison between method-resource combinations

1. Run `*_prep.R` scripts to generate Seurat objects.   

2. Run `main_bash.sh` from bash_scripts folder to generate liana output with all resource and method combinations.  

3. Run `comparison_analyse.R` to generate data for the manuscript.  

4. Run `comparison_plot.R` to generate the plots for the manuscript.  


## Datasets
#### CRC
The processed single cell RNA-Seq data 23 for 23 Korean colorectal cancer patients were obtained from
[GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465)

#### CBMC and panc8
The labelled scRNA-Seq data for pancreatic islet 56 and cord blood mononuclear cells 57 were obtained via SeuratData,
normalized with Seurat, and used for CCC inference without any further formatting and filtering.

#### BRCA
The processed and annotated Human Breast Cancer single-cell atlas (Wu et al. 2021) was obtained from [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078). 
