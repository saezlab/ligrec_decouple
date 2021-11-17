# Script to Prepare Panc8 data for LIANA
require(tidyverse)
require(magrittr)
require(liana)

# Get data from SeuratData
SeuratData::InstallData("panc8")

# Load panc8 Data
panc8data <- SeuratData::LoadData("panc8") %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures()

# Format metadata
panc8data@meta.data %<>%
    mutate(celltype = as.factor(celltype))
panc8data <- subset(panc8data, cells = rownames(panc8data@meta.data))
panc8data <- Seurat::SetIdent(panc8data, value = panc8data@meta.data$celltype)
Seurat::Idents(panc8data) <- panc8data@meta.data$celltype

saveRDS(panc8data, "data/input/comparison/panc8.RDS")
