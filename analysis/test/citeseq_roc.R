# Load prerequisites
require(yardstick)
require(SingleCellExperiment)
require(Seurat)
require(SeuratDisk)
require(liana)
require(tidyverse)
require(magrittr)
require(ComplexHeatmap)
source("analysis/citeseq/citeseq_src.R")


op_resource <- select_resource("OmniPath")[[1]]
murine_resource <- select_resource("OmniPath")[[1]] %>%
    convert_to_murine()
arbitrary_thresh = 1.645 # one-tailed alpha = 0.05
citeseq_dir <- "data/input/citeseq/"


# test ----
seurat_object <- readRDS("data/input/citeseq/10k_pbmcs/10k_pbmcs_seurat.RDS")
liana_res <- readRDS("data/input/citeseq/10k_pbmcs/10k_pbmcs-liana_res-0.1.RDS")
cluster_key = "seurat_clusters"

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
    prepare_for_roc(arbitrary_thresh = arbitrary_thresh)


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

pr_roc_tibble <- readRDS("data/output/citeseq_out/citeseq_aurocs.RDS")

get_auroc_heat(pr_roc_tibble, "prc",
               heatmap_legend_param = list(title="PRAUC"))

get_auroc_heat(pr_roc_tibble, "roc",
               heatmap_legend_param = list(title="AUROC"))





#' ADT-Specificity ROC Helper/Loader Function
#' @inheritParams generate_specificity_roc
#'
#' @returns a tibble with summarized ADT-LR correlations
#'
#' @details Simply loads the needed seurat and liana objects and calls the
#' `wrap_adt_corr` function/pipeline.
#' `seurat_object_path` = seurat_object.RDS path
#' `liana_res_path` = liana_res.RDS path
load_adt_roc <- function(subdir = subdir,
                         dir,
                         sobj_pattern = "_seurat.RDS",
                         liana_pattern,
                         op_resource,
                         cluster_key,
                         organism,
                         arbitrary_thresh,
                         ...){

    message(str_glue("Calculating ROC summary for {subdir}"))

    seurat_object_path <- list.subfiles(subdir = subdir,
                                        dir = citeseq_dir,
                                        pattern = sobj_pattern)


    liana_res_path <- list.subfiles(subdir = subdir,
                                    dir = citeseq_dir,
                                    pattern = liana_pattern)

    generate_specificity_roc(seurat_object = readRDS(seurat_object_path),
                             liana_res = readRDS(liana_res_path),
                             op_resource = op_resource,
                             cluster_key = cluster_key,
                             organism=organism,
                             arbitrary_thresh = arbitrary_thresh,
                             ...)
}




#' Helper function used to prepare `adt_rank` elements
#'
#' @param df `get_rank_adt()` output.
#' @param arbitrary_thresh z-score threshold to calculate ROCs
#'
#' @return tidy data frame with meta information for each experiment and the
#'   response and the predictor value which are required for ROC and
#'   PR curve analysis
prepare_for_roc = function(df, arbitrary_thresh){
    df %>%
        filter(!(method_name %in% c("mean_rank", "median_rank"))) %>%
        # WE NOW REVERT RANKS -> TO BE CHANGED TO SCORES
        dplyr::rename(predictor = value) %>%
        mutate(predictor = predictor * -1) %>%
        group_by(method_name) %>%
        dplyr::mutate(response = case_when(adt_scale > arbitrary_thresh ~ 1,
                                           adt_scale <= arbitrary_thresh ~ 0)) %>%
        mutate(response = factor(response, levels = c(1, 0))) %>%
        dplyr::select(source.target.entity, method_name,
                      adt_scale, predictor, response)
}









#' Function to get Ranks for each method and ADT values
#' to reformat the `wrap_adt_corr` function
get_rank_adt <- function(seurat_object,
                         liana_res,
                         op_resource,
                         cluster_key,
                         organism){

    # Get Symbols of Receptors from OP
    receptor_syms <- str_to_symbol(op_resource$target_genesymbol, organism)

    # convert to singlecell object
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays=list(counts = GetAssayData(seurat_object, assay = "ADT", slot = "counts"),
                    data = GetAssayData(seurat_object, assay = "ADT", slot = "data")),
        colData=DataFrame(label=seurat_object@meta.data[[cluster_key]])
    )
    # Symbols from ADT assay
    adt_symbols <- rownames(sce)

    # Obtain adt genesymbols
    if(organism == "mouse"){
        # obtain both mouse and human but convert to title
        adt_aliases <- get_adt_aliases(adt_symbols = adt_symbols,
                                       organism = "mouse") %>%
            bind_rows(get_adt_aliases(adt_symbols = adt_symbols,
                                      organism = "human")) %>%
            mutate(across(everything(), str_to_title)) %>%
            distinct()
    } else{
        adt_aliases <- get_adt_aliases(adt_symbols = adt_symbols,
                                       organism = organism)
    }

    # check if all matched
    mismatched_adts <- adt_symbols[!(str_to_symbol(adt_symbols, organism) %in%
                                         str_to_symbol(adt_aliases$adt_symbol, organism))]
    message("Missing (ideally none): ",
            glue::glue_collapse(mismatched_adts, sep = ", "))



    # Get ADT stats and bind to LIANA res
    message("Calculating ADT Means")
    adt_means <- get_adt_means(sce = sce,
                               receptor_syms = receptor_syms,
                               adt_aliases = adt_aliases,
                               organism = organism)

    # liana_format_adt
    liana_adt <- liana_res %>%
        filter(receptor %in% adt_means$entity_symbol) %>%
        dplyr::select(source, ligand, target, receptor, ends_with("rank")) %>%
        pivot_longer(c(ligand,receptor),
                     names_to = "type",
                     values_to = "entity_symbol") %>%
        filter(type == "receptor") %>%
        distinct() %>%
        left_join(adt_means, by = c("target", "entity_symbol")) %>% #%>%
        # dplyr::select(ends_with("rank"), starts_with("adt")) %>%
        # pivot_longer(-c(adt_scale, adt_mean),
        #              names_to = "method_name").
        unite(col=source.target.entity, source, target, entity_symbol, sep = ".")  %>%
        dplyr::select(ends_with("rank"), starts_with("adt"), source.target.entity) %>%
        pivot_longer(-c(source.target.entity, adt_scale, adt_mean),
                     names_to = "method_name")

    return(liana_adt)
}




#' Generate ROC
generate_specificity_roc <- function(seurat_object,
                                     liana_res,
                                     op_resource,
                                     cluster_key,
                                     organism,
                                     ...){
    # get adt_ranks DF (need to change it to ADT-SCORE)
    adt_rank_results <- get_rank_adt(seurat_object = seurat_object,
                                     liana_res = liana_res,
                                     op_resource = op_resource,
                                     cluster_key = cluster_key,
                                     organism = organism)

    # To appropriate format (prepare_for_roc) and do ROC
    adt_rank_roc <- adt_rank_results %>%
        prepare_for_roc(arbitrary_thresh = arbitrary_thresh) %>%
        group_by(method_name) %>%
        group_nest(.key = "adt_rank") %>%
        # Calculate ROC
        mutate(roc = .data$adt_rank %>%
                   map(function(df) calc_curve(df))) %>%
        # Calculate PRROC
        mutate(prc = .data$adt_rank %>%
                   map(function(df) calc_curve(df,
                                               curve = "PR",
                                               downsampling = TRUE,
                                               times = 100)))

    return(adt_rank_roc)
}



#' PLOT TO BE FUNCTION -----
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

