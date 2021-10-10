#
require(tidyverse)
require(Seurat)
require(SeuratDisk)

# https://atlas.fredhutch.org/nygc/multimodal-pbmc/
# hfile <- Connect("data/input/seurat_wnn_data/multi.h5seurat")


seurat_object <- LoadH5Seurat("multi.h5seurat",
                              assays="counts", reductions = FALSE,
                              graphs = FALSE, neighbors = FALSE)

saveRDS(seurat_object@meta.data, "metadata.rds")
