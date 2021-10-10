require(tidyverse)
require(magrittr)
require(liana)
require(SPOTlight)
require(Seurat)

# I) Prep Allen Brain Atlas
cortex_sc <- readRDS("data/input/spatial/brain_cortex/allen_cortex.rds")

# Downsample (since SPOTlight shows stable performance with 100 cells,
# we downsample the dataset for computational speed and memory)
cortex_sc@meta.data %<>%
    mutate(subclass = str_replace_all(subclass, "[/]", ".")) %>%
    mutate(subclass = str_replace_all(subclass, " ", ".")) %>%
    rownames_to_column(var = "barcode") %>%
    group_by(subclass) %>%
    slice_sample(n=200) %>%
    ungroup() %>%
    as.data.frame() %>%
    column_to_rownames("barcode")
cortex_sc <- subset(cortex_sc, cells = rownames(cortex_sc@meta.data))
Idents(cortex_sc) <- cortex_sc@meta.data$subclass
gc()

# Normalize appropriately
cortex_sc %<>%
    Seurat::SCTransform(verbose = FALSE) %>%
    Seurat::RunPCA(verbose = FALSE) %>%
    Seurat::RunUMAP(dims = 1:30, verbose = FALSE)

Seurat::DimPlot(cortex_sc,
                group.by = "subclass",
                label = TRUE) + Seurat::NoLegend()

# save the downsampled and formatted object
saveRDS(cortex_sc, "data/input/spatial/brain_cortex/allen_cortex_dwn.rds")





# Deconvolute Mouse slides
c("anterior1", "anterior2", "posterior1", "posterior2")

anterior <- SeuratData::LoadData("stxBrain", type = "anterior1")


anterior2 <- SeuratData::LoadData("stxBrain", type = "anterior2")
posterior <- SeuratData::LoadData("stxBrain", type = "posterior1")
posterior2 <- SeuratData::LoadData("stxBrain", type = "posterior2")





# III) Run LIANA ----
cortex_sc <- readRDS("data/input/spatial/brain_cortex/allen_cortex_dwn.rds")
