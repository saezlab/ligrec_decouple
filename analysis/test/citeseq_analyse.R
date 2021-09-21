require(SingleCellExperiment)
require(Seurat)
require(SeuratDisk)
require(liana)
require(tidyverse)
require(magrittr)

source("analysis/test/citeseq_src.R")


citeseq_dir <- "data/input/citeseq/"
op_resource <- select_resource("OmniPath")[[1]]

# Note that to work this pipeline requires a single {"_seurat.RDS"} object containg
# a v4 seurat_object with RNA and ADT assays
# additionally, it requires {"liana_res"}-containing objects which contain the
# liana results generated for the *single* Seurat object in each subdir

# Iterate over all citeseq directories (i.e. datasets)
list.files(citeseq_dir) %T>%
    # Load .h5 mat and cluster each, and save the appropriate Seurat object
    # map(~load_and_cluster(subdir = .x, dir = citeseq_dir, pattern = ".h5")) %T>%
    # Run liana without expr_prop filtering
    map(~wrap_liana_wrap(subdir = .x, dir = citeseq_dir, expr_prop = 0)) %T>%
    # Run liana with expr_prop filtering
    map(~wrap_liana_wrap(subdir = .x, dir = citeseq_dir, expr_prop = 0.1))


corr_list <- list.files(citeseq_dir) %>%
    map(function(subdir){
        load_adt_lr(dir = citeseq_dir,
                    subdir = subdir,
                    op_resource = op_resource,
                    cluster_key = "seurat_clusters",
                    liana_pattern = "liana_res-0.RDS"
        )
    }) %>% setNames(list.files(citeseq_dir))

corr_list_01 <- list.files(citeseq_dir) %>%
    map(function(subdir){
        load_adt_lr(dir = citeseq_dir,
                    subdir = subdir,
                    op_resource = op_resource,
                    cluster_key = "seurat_clusters",
                    liana_pattern = "liana_res-0.1.RDS"
        )
    }) %>% setNames(list.files(citeseq_dir))

# Plot LR-ADT corr per dataset
corr_list %>%
    enframe() %>%
    unnest(value) %>%
    # remove Mean/Median ranks
    filter(!(method %in% c("median_rank", "mean_rank"))) %>%
    # rename methods
    mutate(method = str_to_title(method)) %>%
    mutate(method = gsub("\\..*","", method)) %>%
    # rename metrics
    mutate(metric = if_else(metric=="scale", "Cluster-specific Mean", metric)) %>%
    mutate(metric = if_else(metric=="mean", "Mean", metric)) %>%
    mutate(metric = if_else(metric=="prop", "Cell Proportion", metric)) %>%
    ggplot(aes(x = metric, y = estimate, colour = method, shape = name)) +
    geom_point(size = 6, position = position_dodge(w = 0.05)) +
    xlab("Expression Type") +
    ylab("Kendal's tau Correlation Coefficients") +
    theme_minimal(base_size = 24) +
    labs(colour = "Method", shape="Dataset")

# plot with expr_prop 0.1
corr_list_01 %>%
    enframe() %>%
    unnest(value) %>%
    # remove Mean/Median ranks
    filter(!(method %in% c("median_rank", "mean_rank"))) %>%
    # rename methods
    mutate(method = str_to_title(method)) %>%
    mutate(method = gsub("\\..*","", method)) %>%
    # rename metrics
    mutate(metric = if_else(metric=="scale", "Cluster-specific Mean", metric)) %>%
    mutate(metric = if_else(metric=="mean", "Mean", metric)) %>%
    mutate(metric = if_else(metric=="prop", "Cell Proportion", metric)) %>%
    ggplot(aes(x = metric, y = estimate, colour = method, shape = name)) +
    geom_point(size = 6, position = position_dodge(w = 0.05)) +
    xlab("Expression Type") +
    ylab("Kendal's tau Correlation Coefficients") +
    theme_minimal(base_size = 24) +
    labs(colour = "Method", shape="Dataset")



