# Script used for the ADT-LR, ADT-RNA correlation (~Proteomics benchmark)

# Note that to work this pipeline requires a single {"_seurat.RDS"} object
# containing a v4 seurat_object with RNA and ADT assays
# additionally, it requires {"liana_res"}-containing objects which contain the
# liana results generated for the *single* Seurat object in each subdir


# Run mouse_citeseq_convert.Rmd to generate seurat objects from the h5anndata
# and format them appropriately for the LIANA analysis
require(liana)
require(tidyverse)
require(magrittr)
require(SingleCellExperiment)
require(Seurat)
require(SeuratDisk)
source("analysis/citeseq/citeseq_src.R")
source("src/eval_utils.R")

citeseq_dir <- "data/input/citeseq/"
op_resource <- select_resource("OmniPath")[[1]]
murine_resource <- readRDS("data/input/murine_omnipath.RDS")


### I) Generate Generate Seurat Objects and run LIANA ----
# Iterate over all citeseq directories (i.e. datasets)
list.files(citeseq_dir) %T>%
    # Load 10x .h5 mat and cluster each, and save the appropriate Seurat object
    map(function(subdir){
        if(subdir %in% c("10k_malt", "10k_pbmcs", "5k_pbmcs", "5k_pbmcs_nextgem")){
            load_and_cluster(subdir = subdir, dir = citeseq_dir, pattern = ".h5")
        } else(
            return()
        )
    }) %T>%
    # Run LIANA
    map(function(subdir){
        # Run LIANA /w mouse specific
        if(stringr::str_detect(subdir, pattern = "spleen")){
            wrap_liana_wrap(subdir = subdir,
                            dir = citeseq_dir,
                            method = c('call_natmi',
                                       'call_connectome', 'logfc', 'cellchat',
                                       'call_sca', 'cellphonedb', "cytotalk"
                            ),
                            squidpy.params=list(cluster_key = "seurat_clusters",
                                                seed = as.integer(1004)),
                            expr_prop = 0.1,
                            cellchat.params = list(nboot=1000,
                                                   expr_prop = 0,
                                                   organism="mouse"),
                            call_natmi.params = list(
                                expr_file = str_glue("{subdir}_em.csv"),
                                meta_file = str_glue("{subdir}_metadata.csv"),
                                output_dir = str_glue("{subdir}_results"),
                                reso_name = str_glue("{subdir}_placeholder")
                            ),
                            resource = "custom",
                            external_resource = murine_resource,
                            organism = "mouse",

                            )
        } else { # human
            wrap_liana_wrap(subdir = subdir,
                            dir = citeseq_dir,
                            method = c('call_natmi',
                                       'call_connectome', 'logfc', 'cellchat',
                                       'call_sca', 'cellphonedb', "cytotalk"
                            ),
                            squidpy.params=list(cluster_key = "seurat_clusters",
                                                seed = as.integer(1004)),
                            expr_prop = 0.1,
                            cellchat.params = list(nboot=1000,
                                                   expr_prop = 0),
                            call_natmi.params = list(
                                expr_file = str_glue("{subdir}_em.csv"),
                                meta_file = str_glue("{subdir}_metadata.csv"),
                                output_dir = str_glue("{subdir}_results"),
                                reso_name = str_glue("{subdir}_placeholder")
                            ),
                            resource = "custom",
                            external_resource = op_resource,
                            organism = "human"
                            )
        }
    })

