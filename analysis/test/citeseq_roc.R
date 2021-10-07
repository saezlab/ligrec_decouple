# Load prerequisites

require(SingleCellExperiment)
require(Seurat)
require(SeuratDisk)
require(liana)
require(tidyverse)
source("analysis/citeseq/citeseq_src.R")


op_resource <- select_resource("OmniPath")[[1]]
murine_resource <- select_resource("OmniPath")[[1]] %>%
    convert_to_murine()
arbitrary_thresh = 1.645 # one-tailed alpha = 0.05
citeseq_dir <- "data/input/citeseq/"


### RUN all sets ----
pr_roc_tibble <- list.files(citeseq_dir) %>%
    # Get ROCS
    map(function(subdir){
        # Run LIANA /w mouse specific
        if(stringr::str_detect(subdir, pattern = "spleen")){
            load_adt_roc(subdir = subdir,
                         dir = citeseq_dir,
                         op_resource = murine_resource,
                         organism = "mouse",
                         cluster_key = "seurat_clusters",
                         arbitrary_thresh = arbitrary_thresh,
                         sobj_pattern = "_seurat.RDS",
                         liana_pattern = "liana_res-0.1.RDS"
                         )
        } else { # human
            load_adt_roc(subdir = subdir,
                         dir = citeseq_dir,
                         op_resource = op_resource,
                         organism = "human",
                         cluster_key = "seurat_clusters",
                         arbitrary_thresh = arbitrary_thresh,
                         sobj_pattern = "_seurat.RDS",
                         liana_pattern = "liana_res-0.1.RDS"
            )
        }
    }) %>% setNames(list.files(citeseq_dir)) %>%
    enframe() %>%
    unnest(value) %>%
    rename(dataset = name)
saveRDS(pr_roc_tibble, "data/output/citeseq_out/citeseq_aurocs.RDS")

# Plots
pr_roc_tibble <- readRDS("data/output/citeseq_out/citeseq_aurocs.RDS")

get_auroc_heat(pr_roc_tibble, "roc",
               heatmap_legend_param = list(title="AUROC"))

get_auroc_heat(pr_roc_tibble, "prc",
               heatmap_legend_param = list(title="PRAUC"))



# test ----
seurat_object <- readRDS("data/input/citeseq/10k_pbmcs/10k_pbmcs_seurat.RDS")
liana_res <- readRDS("data/input/citeseq/10k_pbmcs/10k_pbmcs-liana_res-0.1.RDS")
cluster_key = "seurat_clusters"
organism = "human"

pp_decoupler <- readRDS("~/Downloads/php_result.rds")

test <- generate_specificity_roc(seurat_object = seurat_object,
                                 liana_res = liana_res,
                                 op_resource = op_resource,
                                 cluster_key = "seurat_clusters",
                                 organism="human",
                                 arbitrary_thresh = arbitrary_thresh)


# get adt_ranks DF (need to change it to ADT-SCORE)
adt_rank_results <- get_rank_adt(seurat_object = seurat_object,
                                 liana_res = liana_res,
                                 op_resource = op_resource,
                                 cluster_key = cluster_key,
                                 organism = organism)

# threshold for TP
hist(adt_rank_results$adt_scale)

# To appropriate format (prepare_for_roc)
adt_rank_prep <- adt_rank_results %>%
    prepare_for_roc(arbitrary_thresh = arbitrary_thresh) %>%
    ungroup()

adt_rank_prep %>% calc_curve()


# ROC DF
df_roc <- adt_rank_prep %>%
    # +dataset_name
    dplyr::select(method_name, roc) %>%
    unnest(roc)


# PRROC DF
df_prc <- adt_rank_prep %>%
    # +dataset_name
    dplyr::select(method_name, prc) %>%
    unnest(prc)

# ROC plot
ggplot(df_roc, aes(x = 1-specificity,
                   y = sensitivity,
                   colour = .data$method_name)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")


#' PLOTs raw -----
specificity_summary(adt_rank_roc){
    # ROC DF
    df_roc <- adt_rank_roc %>%
        # +dataset_name
        dplyr::select(method_name, roc) %>%
        unnest(roc)

    # PRROC DF
    df_prc <- adt_rank_roc %>%
        # +dataset_name
        dplyr::select(method_name, prc) %>%
        unnest(prc)

    # ROC plot
    ggplot(df_roc, aes(x = 1-specificity,
                       y = sensitivity,
                       colour = .data$method_name)) +
        geom_line() +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        xlab("FPR (1-specificity)") +
        ylab("TPR (sensitivity)")

    # ROC heatmap
    roc_heat <- df_roc %>%
        as.data.frame() %>%
        dplyr::select(method_name, auc) %>%
        distinct() %>%
        column_to_rownames("method_name") %>%
        pheatmap(.,
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 treeheight_col = 0,
                 treeheight_row = 0,
                 display_numbers = TRUE,
                 silent = TRUE)
    roc_heat

    # PRROC heatmap
    prc_heat <- df_prc %>%
        as.data.frame() %>%
        dplyr::select(method_name, auc) %>%
        distinct() %>%
        column_to_rownames("method_name") %>%
        pheatmap(.,
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 treeheight_col = 0,
                 treeheight_row = 0,
                 display_numbers = TRUE,
                 silent = TRUE)
    prc_heat

}

