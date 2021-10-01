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
    # Load 10x .h5 mat and cluster each, and save the appropriate Seurat object
    # map(~load_and_cluster(subdir = .x, dir = citeseq_dir, pattern = ".h5")) %T>%
    # Run liana without expr_prop filtering
    # map(~wrap_liana_wrap(subdir = .x, dir = citeseq_dir, expr_prop = 0)) %T>%
    # Run liana with expr_prop filtering
    # map(~wrap_liana_wrap(subdir = .x, dir = citeseq_dir, expr_prop = 0.1))
    # Run with the original methods
    # map(~wrap_liana_wrap(subdir = .x,
    #                      dir = citeseq_dir,
    #                      expr_prop = 0.1, # only applied to Squidpy and cellchat
    #                      method = c("call_natmi", "call_connectome", "logfc",
    #                                 "cellchat", "call_sca", "squidpy")))


corr_list <- list.files(citeseq_dir) %>%
    map(function(subdir){
        load_adt_lr(dir = citeseq_dir,
                    subdir = subdir,
                    op_resource = op_resource,
                    cluster_key = "seurat_clusters",
                    liana_pattern = "liana_res-0.RDS"
        )
    }) %>% setNames(list.files(citeseq_dir)) %>%
    enframe() %>%
    unnest(value) %>%
    rename(dataset = name)

corr_list_01 <- list.files(citeseq_dir) %>%
    map(function(subdir){
        # If mouse, load convert use murine-specific conversion
        if(stringr::str_detect(subdir, pattern = "spleen")){
        } else { # human
            load_adt_lr(dir = citeseq_dir,
                        subdir = subdir,
                        op_resource = op_resource,
                        cluster_key = "seurat_clusters",
                        liana_pattern = "liana_res-0.1.RDS"
            )
        }

    }) %>% setNames(list.files(citeseq_dir)) %>%
    enframe() %>%
    unnest(value) %>%
    rename(dataset = name)







# Plot LR-ADT corr per dataset (dotplot)
corr_list %>%
    # remove Mean/Median ranks
    filter(!(method %in% c("median_rank", "mean_rank"))) %>%
    # rename methods
    mutate(method = str_to_title(method)) %>%
    mutate(method = gsub("\\..*","", method)) %>%
    # rename metrics
    mutate(metric = if_else(metric=="scale", "Cluster-specific Mean", metric)) %>%
    mutate(metric = if_else(metric=="mean", "Mean", metric)) %>%
    mutate(metric = if_else(metric=="prop", "Cell Proportion", metric)) %>%
    ggplot(aes(x = metric, y = estimate, colour = method, shape = dataset)) +
    geom_point(size = 6, position = position_dodge(w = 0.05)) +
    xlab("Expression Type") +
    ylab("Kendal's tau Correlation Coefficients") +
    theme_minimal(base_size = 24) +
    labs(colour = "Method", shape="Dataset")



corr_list %>%
    # remove Mean/Median ranks
    filter(!(method %in% c("median_rank", "mean_rank"))) %>%
    # rename methods
    mutate(method = str_to_title(method)) %>%
    mutate(method = gsub("\\..*","", method)) %>%
    # rename metrics
    mutate(metric = if_else(metric=="scale", "Cluster-specific Mean", metric)) %>%
    mutate(metric = if_else(metric=="mean", "Mean", metric)) %>%
    mutate(metric = if_else(metric=="prop", "Cell Proportion", metric)) %>%
    ggplot(aes(x = factor(method), y = estimate)) +
    geom_boxplot(aes(fill = method), alpha = 0.20) +
    geom_jitter(aes(color = method, shape = dataset), size = 4) +
    xlab("") +
    ylab("Kendal's tau Coefficient") +
    theme_minimal(base_size = 24) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(~metric, scales = "free")

corr_list_01 %>%
    # remove Mean/Median ranks
    filter(!(method %in% c("median_rank", "mean_rank"))) %>%
    # rename methods
    mutate(method = str_to_title(method)) %>%
    mutate(method = gsub("\\..*","", method)) %>%
    # rename metrics
    mutate(metric = if_else(metric=="scale", "Cluster-specific Mean", metric)) %>%
    mutate(metric = if_else(metric=="mean", "Mean", metric)) %>%
    mutate(metric = if_else(metric=="prop", "Cell Proportion", metric)) %>%
    ggplot(aes(x = factor(method), y = estimate)) +
    geom_boxplot(aes(fill = method), alpha = 0.20) +
    geom_jitter(aes(color = method, shape = dataset, size = n_genes)) +
    xlab("") +
    ylab("Kendal's tau Coefficient") +
    theme_minimal(base_size = 24) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(~metric, scales = "free")




# CBMC test
seurat_object = readRDS("data/input/citeseq/cmbcs/cbmc_seurat.RDS")
test <- liana_wrap(seurat_object = seurat_object) %>%
    liana_aggregate()
saveRDS(test, "data/input/citeseq/cmbcs/cmbcs-liana_res-0.1.RDS")



### Run on Murine -----
seurat_object <- readRDS("data/input/citeseq/spleen_lymph_206//spleen_lymph_206_seurat.RDS")
op_resource <- liana::select_resource("OmniPath")[[1]] %>%
    convert_to_murine()

murine_liana_206 <- liana_wrap(seurat_object,
                               resource = "custom",
                               external_resource = op_resource,
                               expr_prop = 0.1,
                               cellchat.params=list(organism="mouse"))

# everything to title (some that are nottotitle will be mismatched otherwise, due to the upper of squidpy)
murine_liana %<>% map(function(res) res %>% mutate_at(.vars = c("ligand", "receptor"), str_to_title))

murine_liana %<>% liana_aggregate
saveRDS(murine_liana_206, "data/input/citeseq/spleen_lymph_206/spleen_lymph_206-liana_res-0.1.RDS")



# 111 ----
seurat_object <- readRDS("data/input/citeseq/spleen_lymph_101/spleen_lymph_111_seurat.RDS")
op_resource <- liana::select_resource("OmniPath")[[1]] %>%
    convert_to_murine()

murine_liana_111 <- liana_wrap(seurat_object,
                               resource = "custom",
                               external_resource = op_resource,
                               expr_prop = 0.1,
                               cellchat.params=list(organism="mouse"))

murine_liana_111$squidpy %<>%
    mutate_at(.vars = c("ligand", "receptor"), str_to_title)
murine_liana_111 %<>% liana_aggregate
saveRDS(murine_liana_111, "data/input/citeseq/spleen_lymph_101//spleen_lymph_111-liana_res-0.1.RDS")

### Murine Corr
seurat_object <- readRDS("data/input/citeseq/spleen_lymph_206/spleen_lymph_206_seurat.RDS")
liana_res <- readRDS("data/input/citeseq/spleen_lymph_206/spleen_lymph_206-liana_res-0.1.RDS")
murine_resource <- liana::select_resource("OmniPath")[[1]] %>%
    convert_to_murine()
cluster_key = "seurat_clusters"
organism = "mouse"

corrx <- wrap_adt_corr(seurat_object = seurat_object,
                    liana_res = liana_res,
                    op_resource = murine_resource,
                    cluster_key = "seurat_clusters", # constently using to this name
                    organism = "mouse"
                    )

seurat_object <- readRDS("data/input/citeseq/spleen_lymph_101/spleen_lymph_111_seurat.RDS")
liana_res <- readRDS("data/input/citeseq/spleen_lymph_101/spleen_lymph_111-liana_res-0.1.RDS")
murine_resource <- liana::select_resource("OmniPath")[[1]] %>%
    convert_to_murine()
cluster_key = "seurat_clusters"
organism = "mouse"

corry <- wrap_adt_corr(seurat_object = seurat_object,
                       liana_res = liana_res,
                       op_resource = murine_resource,
                       cluster_key = "seurat_clusters", # constently using to this name
                       organism = "mouse"
                       )
