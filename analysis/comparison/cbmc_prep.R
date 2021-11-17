# Script to Prepare Panc8 data for LIANA
require(tidyverse)
require(magrittr)
require(liana)

# Load from SeuratData
cbmcdata <- SeuratData::LoadData("cbmc")

# Format Metadata
cbmcdata@meta.data %<>%
    mutate(celltype = as.factor(rna_annotations)) %>%
    # remove + as it breaks Squidpy
    mutate(rna_annotations = str_replace_all(rna_annotations, "[+]", "")) %>%
    mutate(rna_annotations = str_replace(rna_annotations, "/", ".")) %>%
    mutate(rna_annotations = if_else(rna_annotations=="B", "B cell", rna_annotations)) %>%
    # filter multiplets and mouse droplets
    filter(!(rna_annotations %in% c("Mouse", "Multiplets", "T Mono doublets"))) %>%
    mutate(rna_annotations = str_replace(rna_annotations, " ", ".")) %>%
    # name to seurat_cluster for consistency with other cmbc seurat objects
    mutate(seurat_clusters = as.factor(rna_annotations))

cbmcdata <- subset(cbmcdata, cells = rownames(cbmcdata@meta.data))
cbmcdata <- Seurat::SetIdent(cbmcdata, value = cbmcdata@meta.data$seurat_clusters)
Seurat::Idents(cbmcdata) <- cbmcdata@meta.data$seurat_clusters

# Basic preprocessing
cbmcdata <- cbmcdata %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures()

# Format Meta
saveRDS(cbmcdata, "data/input/citeseq/cmbcs/cbmc_seurat.RDS")


