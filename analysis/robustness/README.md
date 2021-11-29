# Analyzing Robustness of CCI Inference Methods

## Objectives:

This project aims to assess reliably cell-cell-communication/interaction inference methods based on single cell transcriptomic data. One approach dilutes the resource used by the methods with unvetted ligand-receptor pairs that create false positives (or more accurately non-canonical positives). This imitates both a poorly researched resource, as well as the difference in predictions when switching between resources. The other approach lies in manipulating the cluster annotations of the data set used. By reshuffling or sub-setting a proportion of the cell clusters, we mimic poorly cluster data, or data with a low number of cells.

In both cases, we establish a baseline of highly ranked CCI predictions in normal circumstances and then compare it to highly ranked CCI predictions in the manipulated (diluted, reshuffled, or subset) circumstance. We track and plot the overlap of this pairwise comparison across multiple degrees of manipulation.

## Prerequisites:

This code has all the same dependencies as LIANA (mainly [Seurat](https://satijalab.org/seurat/), [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html), [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment), [base tidyverse](https://www.tidyverse.org/packages/) ) and requires LIANA++ to run. In addition, it makes use of [lubridate](https://lubridate.tidyverse.org/) for runtime calculations. You can get an overview of LIANA [here](https://github.com/saezlab/liana), and learn more about LIANA++ [here](https://saezlab.github.io/liana/articles/liana_devel.html).

## Sub-folders:

1.  **Code**: All the code can be found here, including the the GetData script.
2.  **Data**: When you run the GetData script input data will be collected and deposited here.
3.  **Outputs**: This is where scripts will drop their results if they're saving them to the PC.
