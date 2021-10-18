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


# Run with major celltypes
map2(names(visium_dict),
     visium_dict,
     ~deconv_brca_slides(slide_name = .x,
                         slide_subtype =.y,
                         cluster_key = "celltype_major",
                         n_cells = 100))

# Run with minor cell types
map2(names(visium_dict),
     visium_dict,
     ~deconv_brca_slides(slide_name=.x,
                         slide_subtype =.y,
                         cluster_key = "celltype_minor",
                         n_cells = 100))
