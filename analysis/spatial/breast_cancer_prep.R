# Author: Daniel Dimitrov
# Date: 14.10.2021
# Purpose: Here, I'm preprocessing the Breast cancer (BRCA) single-cell atlas from
# Wu et al., 2021 https://doi.org/10.1038/s41588-021-00911-1, along with the
# visium BRCA slides used in the same publication.

# Load libs
require(tidyverse)
require(magrittr)
require(Seurat)
require(liana)
require(SPOTlight)

# Breast Cancer analysis directory
brca_dir <- "data/input/spatial/Wu_etal_2021_BRCA"


### I) Preprocess scAtlas data to Seurat Object(s) ----
# Read Atlas Expression Matrix
emat <-
    ReadMtx(mtx = file.path(brca_dir, "sc_atlas", "count_matrix_sparse.mtx"),
            cells = file.path(brca_dir, "sc_atlas", "count_matrix_barcodes.tsv"),
            features = file.path(brca_dir, "sc_atlas","count_matrix_genes.tsv"),
            feature.column = 1)

# Create Object
seurat_object <- Seurat::CreateSeuratObject(emat)

# Read and Add Metadata
metadata <- read.csv(file.path(brca_dir, "sc_atlas","metadata.csv"), header = TRUE, row.names = 1)
seurat_object <- Seurat::AddMetaData(seurat_object, metadata)

# Format metadata
seurat_object <- readRDS(file.path(brca_dir, "brca_seurat_atlas.RDS"))
seurat_object@meta.data <- seurat_object@meta.data %>%
    mutate(across(c(celltype_major, celltype_minor, celltype_subset, subtype), ~str_replace_all(.x,"[/]", "\\."))) %>%
    mutate(across(c(celltype_major, celltype_minor, celltype_subset, subtype), ~str_replace_all(.x," ", "\\."))) %>%
    mutate(across(c(celltype_major, celltype_minor, celltype_subset, subtype), ~str_replace_all(.x, "[+]", ""))) %>%
    mutate(across(c(celltype_major, celltype_minor, celltype_subset, subtype), ~str_replace_all(.x, "-", "\\."))) %>%
    mutate(across(c(celltype_major, celltype_minor, celltype_subset, subtype), ~str_replace_all(.x, "_", "\\."))) %>%
    mutate_if(is.character, as.factor)
# Save whole Seurat object
saveRDS(seurat_object, file.path(brca_dir, "brca_atlas_seurat.RDS"))

## Split atlas by BRCA clinical subtype for deconvolution and LIANA =====
seurat_object <- readRDS(file.path(brca_dir, "brca_atlas_seurat.RDS"))
# check if ok
seurat_object@meta.data %>% as_tibble()

# split by subtype
subtypes <- seurat_object@meta.data$subtype %>% levels()
map(subtypes, function(subtype){
    # metadata filtered by clinical subtype
    submeta <- seurat_object@meta.data %>%
        filter(subtype %in% !!subtype)

    # Filter seurat object to clinical subtype cells
    sobj <- subset(seurat_object, cells = rownames(submeta))
    # reassign metadata
    sobj@meta.data <- submeta

    # path to save
    sub_path <- file.path(brca_dir, str_glue("brca_{subtype}_seurat.RDS"))
    # save obj
    saveRDS(sobj, sub_path)

    gc()

    return(submeta %>% as_tibble())
})
