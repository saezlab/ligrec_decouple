require(tidyverse)
require(liana)
require(Seurat)

# Get Args from std in
args <- commandArgs(trailingOnly=TRUE)

# Get Job name
job_name <- args[[2]] # e.g. crc
seurat_path <- args[[3]] # e.g. data/input/comparison/crc_data/crc_korean_form.rds

# Script name
message(str_glue("Running: {args[[1]]} for {job_name}"))

# All resources - Default
resources <- select_resource("all")[-1] %>% names()
message(str_glue("Resources: {glue::glue_collapse(resources, sep = '_')}"))

# Load Seurat_Object
message(str_glue("Reading: {seurat_path}"))
seurat_object <- readRDS(seurat_path)
message(str_glue("Default Assay: {Seurat::DefaultAssay(seurat_object)}"))

# Here, we run LIANA but with multiple resources
liana_res <- liana_wrap(seurat_object,
                        method = c('call_natmi',
                                   'call_connectome', 'logfc', 'cellchat',
                                   'call_sca', 'cellphonedb', "cytotalk"),
                        resource = resources,
                        # this is passed only to squidpy, cellchat, cytotalk, and logfc
                        expr_prop=0.1,
                        cellchat.params = list(nboot=1000,
                                               expr_prop = 0),
                        call_natmi.params = list(
                            expr_file = str_glue("{job_name}_em.csv"),
                            meta_file = str_glue("{job_name}_metadata.csv"),
                            output_dir = str_glue("{job_name}_results"),
                            reso_name = str_glue("{job_name}_placeholder"),
                            .overwrite_data = FALSE),
                        assay = Seurat::DefaultAssay(seurat_object))

# save LIANA results
message("Saving LIANA RESULTS")
saveRDS(liana_res,
        file.path("data/output/comparison_out/",
                  str_glue("{job_name}_liana_res.RDS")))
