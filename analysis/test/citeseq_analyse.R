# Script used for the ADT-LR, ADT-RNA correlation (~Proteomics benchmark)

# Note that to work this pipeline requires a single {"_seurat.RDS"} object
# containing a v4 seurat_object with RNA and ADT assays
# additionally, it requires {"liana_res"}-containing objects which contain the
# liana results generated for the *single* Seurat object in each subdir


# Run mouse_citeseq_convert.Rmd to generate seurat objects from the h5anndata
# and format them appropriately for the LIANA analysis
require(SingleCellExperiment)
require(Seurat)
require(SeuratDisk)
require(liana)
require(tidyverse)
require(magrittr)
source("analysis/test/citeseq_src.R")

citeseq_dir <- "data/input/citeseq/"
op_resource <- select_resource("OmniPath")[[1]]
murine_resource <- op_resource %>%
    convert_to_murine()

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
                            expr_prop = 0.1, # only applied to Squidpy and cellchat
                            # method = c("call_natmi", "call_connectome", "logfc",
                            #            "cellchat", "call_sca", "squidpy"),
                            squidpy.params=list(cluster_key = "seurat_clusters",
                                                seed = as.integer(1)),
                            cellchat.params=list(organism="mouse",
                                                 nboot = 100),
                            resource = "custom",
                            external_resource = murine_resource
                            )
        } else { # human
            wrap_liana_wrap(subdir = subdir,
                            dir = citeseq_dir,
                            expr_prop = 0.1, # only applied to Squidpy and cellchat
                            # method = c("call_natmi", "call_connectome", "logfc",
                            #            "cellchat", "call_sca", "squidpy"),
                            squidpy.params=list(cluster_key = "seurat_clusters",
                                                seed = as.integer(1)),
                            cellchat.params=list(nboot = 100),
                            resource = "custom",
                            external_resource = op_resource)
        }
    }) %>% setNames(list.files(citeseq_dir)) %>%
    enframe() %>%
    unnest(value) %>%
    rename(dataset = name)


corr_table <- list.files(citeseq_dir) %>%
    map(function(subdir){
        # If mouse, load convert use murine-specific conversion
        if(stringr::str_detect(subdir, pattern = "spleen")){
            load_adt_lr(dir = citeseq_dir,
                        subdir = subdir,
                        op_resource = murine_resource,
                        cluster_key = "seurat_clusters",
                        liana_pattern = "liana_res-0.1.RDS",
                        organism = "mouse"
            )
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


corr_table %>%
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
    geom_boxplot(aes(fill = method), alpha = 0.15) +
    geom_jitter(aes(color = method, shape = dataset, size = n_genes))  +
    scale_shape_manual(values = rep(4:12, len = 7)) +
    xlab("") +
    ylab("Kendal's tau Coefficient") +
    theme_minimal(base_size = 24) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(~metric, scales = "free")
