require(SingleCellExperiment)
require(Seurat)
require(SeuratDisk)
require(liana)
require(tidyverse)
require(magrittr)

source("analysis/test/citeseq_src.R")


citeseq_dir <- "data/input/citeseq/"

# Iterate over all citeseq directories (i.e. datasets)
list.files(citeseq_dir) %T>%
    # Load .h5 mat and cluster each, and save the appropriate Seurat object
    map(~load_and_cluster(subdir = .x, dir = citeseq_dir)) %T>%
    # Run liana without expr_prop filtering
    map(~wrap_liana_wrap(subdir = .x, dir = citeseq_dir, expr_prop = 0)) %T>%
    # Run liana with expr_prop filtering
    map(~wrap_liana_wrap(subdir = .x, dir = citeseq_dir, expr_prop = 0.1))

list.files(citeseq_dir) %>%
    map(function(subdir){
        list.subfiles(subdir = subdir,
                      dir = citeseq_dir,
                      pattern = "liana_res")
    }) %>%
    flatten()





# Load files
get_all_res <- list.files(citeseq_dir) %>%
    map(function(subdir){
        list.subfiles(subdir = subdir,
                      dir = citeseq_dir,
                      pattern = "liana_res")
    }) %>%
    flatten() %>%
    map(function(liana_obj_path){
        # liana_res file name with dataset and proportion
        name_of_file <- str_split(liana_obj_path, pattern = "/")[[1]][6] %>%
            substr(., 1, nchar(.)-4) %>%
            as.character()

        list(name_of_file = readRDS(liana_obj_path)) %>%
            setNames(name_of_file)
    }) %>%
    flatten()

#



