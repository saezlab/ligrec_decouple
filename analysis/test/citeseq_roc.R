# Load prerequisites
require(yardstick)
require(SingleCellExperiment)
require(Seurat)
require(SeuratDisk)
require(liana)
require(tidyverse)
require(magrittr)
require(pheatmap)
source("analysis/citeseq/citeseq_src.R")


# test
seurat_object <- readRDS("data/input/citeseq/10k_pbmcs/10k_pbmcs_seurat.RDS")
liana_res <- readRDS("data/input/citeseq/10k_pbmcs/10k_pbmcs-liana_res-0.1.RDS")
op_resource <- select_resource("OmniPath")[[1]]# %>%
    # convert_to_murine()
cluster_key = "seurat_clusters"
organism="human"
pp_decoupler <- readRDS("~/Downloads/php_result.rds")
arbitrary_thresh = 1.282 # alpha = 0.1



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





#' Helper function used to prepare `adt_rank` elements
#'
#' @param df `activity` column elements - i.e. `get_rank_adt()` output.
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
                      adt_scale, predictor, response) %>%
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
}




#' Helper function to produce AUROC heatmap
#' @param auroc_tibble Tibble with calculated AUROC
#' @return returns an AUROC or Precision-Recall AUC heatmap
#' @import pheatmap ggplot2
get_auroc_heat <- function(auroc_tibble){

}





#' This function takes the elements of the `activity` column and calculates
#'    precision-recall and ROC curves (depending on `curve`).
#' The `activity` column is populated with the output for each stat method and
#'    results from the `run_benchmark()` function. Each of the elements
#'    in `activity` are results from runs of the \link{decouple} wrapper.
#'
#' @param df run_benchmark roc column provided as input
#' @param downsampling logical flag indicating if the number of Negatives
#'    should be downsampled to the number of Positives
#' @param times integer showing the number of downsampling
#' @param curve whether to return a Precision-Recall Curve ("PR") or ROC ("ROC")
#' @param seed An integer to set the RNG state for random number generation. Use
#'    NULL for random number generation.
#'
#' @return tidy data frame with precision, recall, auc, n, cp, cn and coverage
#'    in the case of PR curve; or sensitivity and specificity, auc, n, cp, cn
#'    and coverage in the case of ROC.
#' @import yardstick
#'
#' @export
calc_curve = function(df,
                      downsampling = FALSE,
                      times = 1000,
                      curve = "ROC",
                      seed = 420){
    set.seed(seed)

    if(curve=="PR"){
        res_col_1 <- "precision"
        res_col_2 <- "recall"
        curve_fun = yardstick::pr_curve
        auc_fun = yardstick::pr_auc
    }
    else if(curve=="ROC"){
        res_col_1 <- "sensitivity"
        res_col_2 <- "specificity"
        curve_fun = yardstick::roc_curve
        auc_fun = yardstick::roc_auc
    }

    if (sum(which(df$response == 0)) == nrow(df)){
        return(as_tibble(NULL))
    }

    cn = df %>% filter(.data$response == 0)
    cp = df %>% filter(.data$response == 1)

    feature_coverage = length(unique(df$source.target.entity))

    if(downsampling == TRUE){
        num_tp = nrow(cp)

        res = map_df(seq(from=1, to=times, by=1), function(i) {
            df_sub = sample_n(cn, num_tp, replace=TRUE) %>%
                bind_rows(cp)

            r_sub = df_sub %>%
                curve_fun(.data$response, .data$predictor)

            auc = df_sub %>%
                auc_fun(.data$response, .data$predictor) %>%
                pull(.data$.estimate)

            res_sub = tibble({{ res_col_1 }} := r_sub %>% pull(res_col_1),
                             {{ res_col_2 }} := r_sub %>% pull(res_col_2),
                             th = r_sub$.threshold,
                             auc = auc,
                             n = length(which(df$response == 1)),
                             cp = nrow(cp),
                             cn = nrow(cn),
                             coverage = feature_coverage) %>%
                mutate("run" = i)

        })
        # Get Average AUC
        res <- res %>% dplyr::rename("raw_auc" = auc)
        # auc is the mean of all iterations, raw_auc is the value per iteration
        res$auc <- sum(res$raw_auc)/length(res$raw_auc)
        res$cn <- nrow(cp)

    } else {
        r = df %>%
            curve_fun(.data$response, .data$predictor)
        auc = df %>%
            auc_fun(.data$response, .data$predictor)

        res = tibble({{ res_col_1 }} := r %>% pull(res_col_1),
                     {{ res_col_2 }} := r %>% pull(res_col_2),
                     th = r$.threshold,
                     auc = auc$.estimate,
                     n = length(which(df$response == 1)),
                     cp = nrow(cp),
                     cn = nrow(cn),
                     coverage = feature_coverage) %>%
            arrange(!!res_col_1, !!res_col_2)
    }
    return(res)
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

