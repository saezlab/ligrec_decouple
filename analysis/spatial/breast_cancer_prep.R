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
# Run function to subsample atlases accordingly
map(c("ER", "TNBC"),
    function(slide_subtype,
             cluster_key = "celltype_minor",
             n_cells = 100){ # 500 for major, 200 for minor


        # Define deconvolution results and input dir
         deconv_directory <- file.path(project_dir, brca_dir,
                                       "deconv", str_glue("{slide_subtype}_{cluster_key}"))
         dir.create(deconv_directory, showWarnings = FALSE)
         message(str_glue("Created: {deconv_directory}"))

         # Load BRCA Subtype Atlas
         atlas_object <- readRDS(file.path(project_dir,
                                           brca_dir,
                                           str_glue("brca_{slide_subtype}_seurat.RDS")))

         # subsample
         message("Subsampling")
         set.seed(1234)
         submeta <- atlas_object@meta.data %>%
             rownames_to_column("barcode") %>%
             group_by(!!sym(cluster_key)) %>%
             mutate(ncels_by_group = n()) %>%
             filter(ncels_by_group >= 25) %>% # as in Wu et al., 2021
             slice_sample(n=n_cells) %>% # downsample, SPOTlight stable with 100 cells
             #mutate(!!cluster_key = as.factor(as.character(cluster_key)))
             ungroup() %>%
             mutate({{ cluster_key }} := as.factor(as.character(.data[[cluster_key]]))) %>%
             as.data.frame() %>%
             column_to_rownames("barcode")
         # reassign meta/clusters
         atlas_object@meta.data <- submeta
         Idents(atlas_object) <- submeta[[cluster_key]]
         seurat_object[[cluster_key]] <- Idents(seurat_object)
         # subset to meta
         atlas_object <- subset(atlas_object, cells = rownames(submeta))

         # Normalize object
         atlas_object %<>%
             Seurat::SCTransform(verbose = FALSE) %>%
             Seurat::RunPCA(verbose = FALSE) %>%
             Seurat::RunUMAP(dims = 1:30, verbose = FALSE)
         saveRDS(atlas_object, file.path(deconv_directory,
                                         str_glue("{slide_subtype}_{cluster_key}_sub_seurat.RDS")))

         # Find Markers
         # cluster_markers_all <-
         #     Seurat::FindAllMarkers(object = atlas_object,
         #                            assay = "SCT",
         #                            slot = "data",
         #                            verbose = TRUE,
         #                            only.pos = TRUE)
         # saveRDS(cluster_markers_all,
         #         file.path(deconv_directory,
         #                   str_glue("{slide_subtype}_{cluster_key}_markers.RDS")))
         return()
})

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
deconv_spec_plots$`ER_celltype_minor/CID4290_ER_celltype_minor_deconv.RDS`[[2]]
deconv_spec_plots$`ER_celltype_major/CID4290_ER_celltype_major_deconv.RDS`[[2]]



### RUN LIANA on ER and TNBC subtype atlases ----
# Tibble with all subtype - celltype combinations
all_combs <- tibble(slide_subtype = c("TNBC", "ER", "TNBC", "ER"),
                    cluster_key = c("celltype_major", "celltype_minor",
                                    "celltype_minor", "celltype_major"))

pmap(all_combs, function(slide_subtype, cluster_key){
    deconv_directory <- file.path(brca_dir,
                                  "deconv",
                                  str_glue("{slide_subtype}_{cluster_key}"))
    message(deconv_directory)

    seurat_object <- readRDS(file.path(deconv_directory,
                                       str_glue("{slide_subtype}_{cluster_key}_sub_seurat.RDS")))

    # RUN LIANA
    liana_res <- liana_wrap(seurat_object,
                            assay = "SCT",
                            expr_prop=0.1,
                            squidpy.params = list(cluster_key = cluster_key))
    saveRDS(liana_res,
            file.path(deconv_directory,
                      str_glue("{slide_subtype}_{cluster_key}_liana_res.RDS")))
})
