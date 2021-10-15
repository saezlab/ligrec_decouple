require(tidyverse)
require(magrittr)
require(Seurat)
require(liana)
require(SPOTlight)
source("analysis/spatial/spatial_src.R")


# Breast Cancer analysis directory
brca_dir <- "data/input/spatial/Wu_etal_2021_BRCA"
project_dir <- getwd()

visium_dict <- list("1142243F" = "TNBC",
                    "1160920F" = "TNBC",
                    "CID4290" = "ER",
                    "CID4465" = "TNBC",
                    "CID4535" = "ER",
                    "CID44971" = "TNBC")

# Run major cell types with 500 cells per cluster -----
map2(names(visium_dict), visium_dict,
     .f=function(slide_name,
                 slide_sybtype,
                 cluster_key = "celltype_minor",
                 n_cells = 100){

         # deconvolution subdirectory
         deconv_directory <- file.path(project_dir, brca_dir,
                                       "deconv", str_glue("{slide_sybtype}_{cluster_key}"))
         message(str_glue("Now running {slide_name}({slide_sybtype}) using {cluster_key} with {n_cells} cells"))
         message(str_glue("To be Saved/Loaded in: {deconv_directory} "))

         # Load Subsampled atlas
         atlas_object <-
             readRDS(file.path(deconv_directory,
                               str_glue("{slide_sybtype}_{cluster_key}_sub_seurat.RDS")))

         # marker genes for that atlas
         cluster_markers_all <-
             readRDS(file.path(deconv_directory,
                       str_glue("{slide_sybtype}_{cluster_key}_markers.RDS")))

         # load visium object
         vis_obj <- Load10X_Spatial_enh(file.path(project_dir,
                                                  brca_dir,
                                                  "/visium_rearranged",
                                                  slide_name),
                                        project_dir = project_dir)

         # Run deconvolution as from tutorial and paper
         set.seed(1234)
         spotlight_ls <- spotlight_deconvolution(
             se_sc = atlas_object,
             counts_spatial = vis_obj@assays$Spatial@counts,
             clust_vr = cluster_key, # cell-type annotations
             cluster_markers = cluster_markers_all, # df with marker genes
             cl_n = n_cells, # number of cells per cell type
             hvg = 3000, # Number of HVG
             ntop = NULL, # How many of the marker genes to use (by default all)
             transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorization and NLS
             method = "nsNMF", # Factorization method
             min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
             )

         saveRDS(spotlight_ls,
                 file.path(deconv_directory,
                            str_glue("{slide_name}_{slide_sybtype}_{cluster_key}_deconv.RDS")))

         return()
})

# # deconvolution subdirectory
# deconv_directory <- file.path(project_dir, brca_dir,
#                               "deconv", str_glue("{slide_sybtype}_{cluster_key}"))
#
# # Load Subsampled atlas
# file.path(deconv_directory,
#           str_glue("{slide_sybtype}_{cluster_key}_sub_seurat.RDS"))
#
# # marker genes for that atlas
# file.path(deconv_directory,
#           str_glue("{slide_sybtype}_{cluster_key}_markers.RDS"))
#
#
# # Load atlas
# seurat_object <- readRDS(file.path(project_dir,
#                                    brca_dir,
#                                    str_glue("brca_{vis_subtype}_seurat.RDS")))
#
# set.seed(1234)
# # subsample
# submeta <- seurat_object@meta.data %>%
#     rownames_to_column("barcode") %>%
#     group_by(!!sym(cluster_key)) %>%
#     mutate(ncels_by_group = n()) %>%
#     filter(ncels_by_group >= 25) %>% # as in Wu et al., 2021
#     slice_sample(n=500) %>% # downsample, SPOTlight stable with 100 cells
#     ungroup() %>%
#     as.data.frame() %>%
#     column_to_rownames("barcode")
# seurat_object <- subset(seurat_object, cells = rownames(submeta))
# Idents(seurat_object) <- seurat_object@meta.data[[cluster_key]]
#
# # Normalize object
# seurat_object %<>%
#     Seurat::SCTransform(verbose = FALSE) %>%
#     Seurat::RunPCA(verbose = FALSE) %>%
#     Seurat::RunUMAP(dims = 1:30, verbose = FALSE)
#
# # Check quickly
# Seurat::DimPlot(seurat_object,
#                 group.by = cluster_key,
#                 label = TRUE) + Seurat::NoLegend()
# saveRDS(seurat_object,
#         file.path(project_dir,
#                   brca_dir,
#                   str_glue("brca_{vis_subtype}_sub_seurat.RDS")))
#
# # Find Markers
# cluster_markers_all <- Seurat::FindAllMarkers(object = seurat_object,
#                                               assay = "SCT",
#                                               slot = "data",
#                                               verbose = TRUE,
#                                               only.pos = TRUE)
# saveRDS(cluster_markers_all,
#         file = file.path(project_dir,
#                          brca_dir,
#                          str_glue("brca_{vis_subtype}_{cluster_key}_markers.RDS")))
#
#
# # Load visium slide
# vis_name <- names(visium_dict)[[1]]
# vis_subtype <- visium_dict[[1]]
# deconv_directory <- file.path(file.path(project_dir, brca_dir, str_glue("deconv_{vis_subtype}_{cluster_key}")))
# deconv_directory
# dir.create(deconv_directory, showWarnings = FALSE)
# names(visium_dict)
#
#
# vis_obj <- Load10X_Spatial_enh(file.path(project_dir,
#                                          brca_dir, "/visium_rearranged",
#                                          vis_name),
#                                project_dir = project_dir)
#
# seurat_object <- readRDS(file.path(project_dir,
#                                    brca_dir,
#                                    str_glue("brca_{vis_subtype}_sub_seurat.RDS")))
# cluster_markers_all <- readRDS(file.path(project_dir,
#                                          brca_dir,
#                                          str_glue("brca_{vis_subtype}_{cluster_key}_markers.RDS")))
#
#
# # Run deconvolution as from tutorial and paper
# set.seed(1234)
# spotlight_ls <- spotlight_deconvolution(
#     se_sc = seurat_object,
#     counts_spatial = vis_obj@assays$Spatial@counts,
#     clust_vr = cluster_key, # cell-type annotations
#     cluster_markers = cluster_markers_all, # df with marker genes
#     cl_n = 200, # number of cells per cell type
#     hvg = 3000, # Number of HVG
#     ntop = NULL, # How many of the marker genes to use (by default all)
#     transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorization and NLS
#     method = "nsNMF", # Factorization method
#     min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
# )

