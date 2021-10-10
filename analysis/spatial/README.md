# Spatial-LR score Evaluation

## Reproduce CiteSeq Analysis

```
```

## Datasets
Spatial transcriptomics datasets (10x visium slides) on sagittal adult mouse brain
anterior and posterior slices were obtained from [`SeuratData`v0.2.1](https://github.com/satijalab/seurat-data),
available under the dataset name of `stxBrain`.

Cell type mapping on the mouse brain slides was performed using [`SPOTlight`](https://academic.oup.com/nar/article/49/9/e50/6129341#248806291), with
100 cells per cluster as in the SPOTlight publication. 

As reference for the cell type mapping, we used [a processed Seurat](https://satijalab.org/seurat/articles/spatial_vignette.html) object of the Allen Brain Atlas ([Tasic et al., 2016](https://www.nature.com/articles/nn.4216))
