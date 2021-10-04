# test run on cluster /w all CRC cell types (needed for correlation /w CRC deconv visium)
require(tidyverse)
require(Seurat)
require(liana)

# makedir



seurat_object <- readRDS("data/input/crc_data/crc_korean.rds")
seurat_object@meta.data <- seurat_object@meta.data %>%
    # Filter Unknown and Unspecified
    filter(!(str_detect("Unknown", Cell_subtype))) %>%
    filter(!(str_starts("Unspecified Plasma", Cell_subtype))) %>%
    # _ breaks Squidpy - need to make it so that it works...
    unite(Cell_type, Cell_subtype, col = "cluster_key", remove = FALSE, sep="-") %>%
    # remove + as it breaks Squidpy
    mutate(cluster_key = str_replace_all(cluster_key, "[+]", "")) %>%
    mutate(cluster_key = str_replace_all(cluster_key, "[/]", "")) %>%
    mutate(cluster_key = factor(cluster_key))# %>%
    # Subset
    # rownames_to_column(var = "barcode") %>%
    # group_by(cluster_key) %>%
    # slice_sample(prop = 0.5) %>%
    # ungroup() %>%
    # as.data.frame() %>%
    # column_to_rownames("barcode")

seurat_object <- subset(seurat_object,
                        cells = rownames(seurat_object@meta.data))
Idents(seurat_object) <- seurat_object@meta.data$cluster_key
gc()

liana_res <- liana_wrap(seurat_object,
                        resource = "OmniPath")
saveRDS(liana_res, "data/output/crc_all_test.RDS")
liana_res <- readRDS("data/output/crc_all_test.RDS")

#
liana_agg <- liana_res %>% liana_aggregate()
