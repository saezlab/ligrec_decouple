#### Run LIANA with all original methods (with defaults args) + multiple resources
# load libs
require(tidyverse)
require(magrittr)
require(Seurat)
require(liana)

# Get Args from std in
args <- commandArgs(trailingOnly=TRUE)

# Get brca_subtpye
brca_subtype <- args[[1]] # "TNBC", "ER", "HER2"
message(str_glue("Now Running: LIANA with {brca_subtype}"))

# All resources - Default
resources <- select_resource("all")[-1] %>% names()
message(str_glue("Resources: {glue::glue_collapse(resources, sep = '_')}"))

# Spatial Deconv Directory (i.e atlas directory)
brca_dir <- "data/input/spatial/Wu_etal_2021_BRCA" # should change to comparison
deconv_directory <- file.path(brca_dir,
                              "deconv", str_glue("{brca_subtype}_celltype_minor"))

# Load BRCA object from Spatial Directory
seurat_object <- readRDS(file.path(deconv_directory,
                                   str_glue("{brca_subtype}_celltype_minor_seurat.RDS")))


# Here, we run LIANA but with multiple resources
liana_res <- liana_wrap(seurat_object,
                        method = c('call_natmi',
                                   'call_connectome', 'logfc', 'cellchat',
                                   'call_sca', 'cellphonedb', "cytotalk"
                                   ),
                        resource = resources,
                        # this is passed only to squidpy, cellchat, cytotalk, and logfc
                        expr_prop=0.1,
                        cellchat.params = list(nboot=1000,
                                               expr_prop = 0),
                        call_natmi.params = list(
                            expr_file = str_glue("{brca_subtype}_em.csv"),
                            meta_file = str_glue("{brca_subtype}_metadata.csv"),
                            output_dir = str_glue("brca_{brca_subtype}_results"),
                            reso_name = str_glue("{brca_subtype}_placeholder"),
                            .overwrite_data = FALSE
                            ),
                        assay = "SCT") # by default as in CellChat

# save LIANA results
message("Saving LIANA RESULTS")
saveRDS(liana_res,
        file.path("data/output/comparison_out/",
                  str_glue("BRCA_{brca_subtype}_liana_res.RDS")))
