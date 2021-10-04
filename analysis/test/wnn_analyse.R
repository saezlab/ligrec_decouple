# WNN
## Check metadata
wnn_meta <- readRDS("data/input/citeseq/seurat_wnn_data/metadata.rds") # Seurat4 is required
unique(wnn_meta$celltype.l1)
unique(wnn_meta$celltype.l2) # 31 - too many. preferably exclude some?
# e.g. rm eryth, Doublet, all naive cell types (i.e. - 3);  merge cDC 1 and 2, Doublet,
unique(wnn_meta$celltype.l3) # way too many

# 228 antibodies (ADTs)
assay <- readRDS("data/input/citeseq/seurat_wnn_data/adt_expression.rds")




sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts = assay),
                                                  colData=DataFrame(label=wnn_meta$celltype.l2))
