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
source("analysis/spatial/spatial_src.R")


# Breast Cancer analysis directory
brca_dir <- "data/input/spatial/Wu_etal_2021_BRCA"
project_dir <- getwd()


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
    message(str_glue("Subsampling atlas to {subtype} subtype"))

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


### II) Prep single-cell atlases -------------
# Run function to format, preprocess, and run markers for atlases accordingly

# Run for Major celltype
map(c("ER", "TNBC"), ~prep_brca_atlases(slide_subtype = .x,
                                        cluster_key = "celltype_major"))

# Run for Minor celltypes
map(c("ER", "TNBC"), ~prep_brca_atlases(slide_subtype = .x,
                                        cluster_key = "celltype_minor"))


#### Run Deconvolution ----
source("analysis/spatial/breast_cancer_deconv.R") # Run deconvolutions

# Check deconv
deconv_dirs <- list.files(path =
                              file.path(brca_dir, "deconv"),
                          recursive = TRUE,
                          pattern = "deconv.RDS")
deconv_dirs

deconv_spec_plots <- deconv_dirs %>%
    map(function(deconv_path){
        deconv_res <- readRDS(file.path(brca_dir, "deconv", deconv_path))
        nmf_mod <- deconv_res[[1]]
        h <- NMF::coef(nmf_mod[[1]])
        rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
        topic_profile_plts <-
            SPOTlight::dot_plot_profiles_fun(
                h = h,
                train_cell_clust = nmf_mod[[2]])
        # Check topics/cell type specificity
        return(topic_profile_plts)
    }) %>%
    setNames(deconv_dirs)


# Check specificity
deconv_spec_plots$`TNBC_celltype_minor/1160920F_TNBC_celltype_minor_deconv.RDS`[[2]]
deconv_spec_plots$`ER_celltype_minor/CID4290_ER_celltype_minor_deconv.RDS`[[2]]



### RUN LIANA on ER and TNBC subtype atlases ----
# Tibble with all subtype - celltype combinations
all_combs <- tibble(slide_subtype = c("TNBC",
                                      "ER",
                                      "TNBC",
                                      "ER"),
                    cluster_key = c("celltype_major",
                                    "celltype_minor",
                                    "celltype_minor",
                                    "celltype_major"))

pmap(all_combs, function(slide_subtype, cluster_key){
    deconv_directory <- file.path(brca_dir,
                                  "deconv",
                                  str_glue("{slide_subtype}_{cluster_key}"))
    message(deconv_directory)

    seurat_object <- readRDS(file.path(deconv_directory,
                                       str_glue("{slide_subtype}_{cluster_key}_seurat.RDS")))
    gc()

    # RUN LIANA
    liana_res <- liana_wrap(seurat_object,
                            assay = "SCT",
                            expr_prop = 0.1,
                            squidpy.params = list(cluster_key = cluster_key))
    saveRDS(liana_res,
            file.path(deconv_directory,
                      str_glue("{slide_subtype}_{cluster_key}_liana_res.RDS")))

    gc()

    return()
})
