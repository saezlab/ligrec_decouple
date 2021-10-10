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
require(yardstick)
require(magrittr)
require(ComplexHeatmap)
source("analysis/citeseq/citeseq_src.R")

citeseq_dir <- "data/input/citeseq/"
op_resource <- select_resource("OmniPath")[[1]]
murine_resource <- readRDS("data/input/murine_omnipath.RDS")
arbitrary_thresh = 1.645 # one-tailed alpha = 0.05


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
                            expr_prop = 0.1, # only applied to Squidpy and cellchat
                            # method = c("call_natmi", "call_connectome", "logfc",
                            #            "cellchat", "call_sca", "squidpy"),
                            squidpy.params=list(cluster_key = "seurat_clusters",
                                                seed = as.integer(1)),
                            cellchat.params=list(organism="mouse",
                                                 nboot = 100),
                            resource = "custom",
                            external_resource = murine_resource,
                            organism = "mouse"
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
                            external_resource = op_resource,
                            organism = "human"
                            )
        }
    }) %>% setNames(list.files(citeseq_dir)) %>%
    enframe() %>%
    unnest(value) %>%
    rename(dataset = name)


### II) Correlations
corr_table <- list.files(citeseq_dir) %>%
    map(function(subdir){
        # If mouse, load convert use murine-specific conversion
        if(stringr::str_detect(subdir, pattern = "spleen")){
            run_adt_pipe(dir = citeseq_dir,
                        subdir = subdir,
                        op_resource = murine_resource,
                        cluster_key = "seurat_clusters",
                        liana_pattern = "liana_res-0.1.RDS",
                        organism = "mouse",
                        adt_pipe_type = "correlation"
            )
        } else { # human
            run_adt_pipe(dir = citeseq_dir,
                        subdir = subdir,
                        op_resource = op_resource,
                        cluster_key = "seurat_clusters",
                        liana_pattern = "liana_res-0.1.RDS",
                        organism = "human",
                        adt_pipe_type = "correlation"
            )
        }

    }) %>% setNames(list.files(citeseq_dir)) %>%
    enframe() %>%
    unnest(value) %>%
    rename(dataset = name)
saveRDS(corr_table, "data/output/citeseq_out/citeseq_correlations.RDS")

# Load and plot
corr_table <- readRDS("data/output/citeseq_out/citeseq_correlations.RDS")
corr_table %<>%
    # remove Mean/Median ranks
    filter(!(method %in% c("median_rank", "mean_rank"))) %>%
    # rename methods
    mutate(method = gsub("\\..*","", method)) %>%
    # recode methods
    mutate(method = recode_methods(method)) %>%
    # recode datasets
    mutate(dataset = recode_datasets(dataset)) %>%
    # rename metrics
    filter(metric=="mean") %>%
    # mutate(metric = if_else(metric=="scale", "Cluster-specific Mean", metric)) %>%
    mutate(metric = if_else(metric=="mean", "Mean", metric)) %>%
    mutate(metric = if_else(metric=="prop", "Cell Proportion", metric))

# Bar plot
corr_table %>%
    ggplot(aes(x = factor(method), y = estimate)) +
    geom_boxplot(aes(fill = method), alpha = 0.15, show.legend = FALSE) +
    geom_jitter(aes(color = method, shape = dataset, size = n_genes))  +
    scale_shape_manual(values = rep(4:12, len = 7)) +
    xlab("") +
    ylab("Kendal's tau Coefficient") +
    theme_minimal(base_size = 24) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(size=guide_legend(title="Receptor Genes per dataset")) +
    labs(fill=guide_legend(title="Method"),
         shape=guide_legend(title="Dataset"))



### III) Receptor Specificity ROC -----
pr_roc_tibble <- list.files(citeseq_dir) %>%
    map(function(subdir){
        # If mouse, load convert use murine-specific conversion
        if(stringr::str_detect(subdir, pattern = "spleen")){
            run_adt_pipe(subdir = subdir,
                         dir = citeseq_dir,
                         op_resource = murine_resource,
                         organism = "mouse",
                         cluster_key = "seurat_clusters",
                         sobj_pattern = "_seurat.RDS",
                         liana_pattern = "liana_res-0.1.RDS",
                         arbitrary_thresh = arbitrary_thresh,
                         adt_pipe_type = "specificity"
            )
        } else { # human
            run_adt_pipe(subdir = subdir,
                         dir = citeseq_dir,
                         op_resource = op_resource,
                         organism = "human",
                         cluster_key = "seurat_clusters",
                         sobj_pattern = "_seurat.RDS",
                         liana_pattern = "liana_res-0.1.RDS",
                         arbitrary_thresh = arbitrary_thresh,
                         adt_pipe_type = "specificity"
            )
        }

    }) %>% setNames(list.files(citeseq_dir)) %>%
    enframe() %>%
    unnest(value) %>%
    rename(dataset = name)

# Save obj
saveRDS(pr_roc_tibble, "data/output/citeseq_out/citeseq_aurocs.RDS")


# Generate Plots
pr_roc_tibble <- readRDS("data/output/citeseq_out/citeseq_aurocs.RDS")

get_auroc_heat(pr_roc_tibble, "roc",
               heatmap_legend_param = list(title="AUROC"))

get_auroc_heat(pr_roc_tibble, "prc",
               heatmap_legend_param = list(title="PRAUC"))
