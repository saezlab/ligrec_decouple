# Script to Prepare CRC data for LIANA
require(tidyverse)
require(magrittr)
require(liana)

#' Helper function to convert CRC data to sparse Seurat
#' @param counts_loc counts location
#' @param meta_loc metadata location
#' @param save_loc location to save sparsified Seurat object
#' @return Return a Sparse Seurat Object
sparsify_to_seurat <- function(counts_loc, meta_loc, save_loc){
    counts <- read_delim(counts_loc,
                         delim = "\t") %>%
        as.data.frame() %>%
        column_to_rownames("Index")

    meta <- read_delim(meta_loc,
                       delim = "\t") %>%
        as.data.frame() %>%
        column_to_rownames("Index")

    crc_seurat <- CreateSeuratObject(counts = Seurat::as.sparse(counts),
                                     project = "10X_CRC") %>%
        Seurat::AddMetaData(meta) %>%
        Seurat::NormalizeData() %>%
        Seurat::FindVariableFeatures() %>%
        saveRDS(., save_loc)
}


# Sparsify CRC data
sparsify_to_seurat(counts_loc = "data/input/crc_data/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt",
                   meta_loc = "data/input/crc_data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt",
                   save_loc = "data/input/crc_data/crc_korean.rds")

# Here we load and format the resulting Seurat object
crc_korean <- readRDS("data/input/comparison/crc_data/crc_korean.rds")

# Format metadata
crc_korean@meta.data %<>%
    # keep only Tumour samples from CMS2
    filter(Class == "Tumor") %>%
    # CMS2 tissues only
    filter(Patient %in% c("SMC07", "SMC09", "SMC11", "SMC18",
                          "SMC21", "SMC22", "SMC23", "SMC25")) %>%
    unite(Cell_type, Cell_subtype, col = "cell_clusters", remove = FALSE, sep=".") %>%
    as.data.frame() %>%
    rownames_to_column("barcode") %>%
    filter(!(str_detect("Unknown", Cell_subtype))) %>%
    mutate(cell_clusters = str_replace_all(cell_clusters, " ", ".")) %>%
    mutate(cell_clusters = gsub("\\+", "", cell_clusters)) %>%
    filter(!(str_starts("Unspecified Plasma", cell_clusters))) %>%
    mutate(cell_clusters = if_else(str_detect(cell_clusters, "SPP1"),
                                   "Myoloids.SPP1+",
                                   cell_clusters)) %>%
    group_by(cell_clusters) %>%
    mutate(cell_num = n()) %>%
    filter(cell_num >= 25) %>%
    ungroup() %>%
    mutate(cell_clusters = as.factor(cell_clusters))  %>%
    as.data.frame(row.names) %>%
    column_to_rownames("barcode")
crc_korean <- subset(crc_korean,
                     cells = rownames(crc_korean@meta.data))
crc_korean <- Seurat::SetIdent(crc_korean,
                               value = crc_korean@meta.data$cell_clusters)

# save formatted Seurat
saveRDS(crc_korean, "data/input/comparison/crc_data/crc_korean_form.rds")

