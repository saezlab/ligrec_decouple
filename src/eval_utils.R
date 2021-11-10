#' Basic function to convert human to mouse genesymbols (temporary solution)
#' @param op_resource omnipath_resource as obtained via `liana::select_resource`
#'
#' @details adapted from https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
convert_to_murine <- function(op_resource){
    require("biomaRt")
    require(liana)
    require(tidyverse)

    # query biomaRt databases
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    # obtain tibble with human and murine genesymbol
    symbols_tibble <- getLDS(attributes = c("hgnc_symbol"),
                             filters = "hgnc_symbol",
                             values = union(op_resource$source_genesymbol,
                                            op_resource$target_genesymbol),
                             mart = human,
                             martL = mouse,
                             attributesL = c("mgi_symbol")) %>%
        dplyr::rename(human_symbol = HGNC.symbol,
                      murine_symbol = MGI.symbol) %>%
        as_tibble()

    # intentionally we introduce duplicates, if needed
    # these should be resolved when LIANA is called
    # as inappropriately matched genes will not be assigned any values
    op_resource %>%
        left_join(symbols_tibble, by=c("target_genesymbol"="human_symbol")) %>%
        mutate(target_genesymbol = murine_symbol, .keep = "unused") %>%
        left_join(symbols_tibble, by=c("source_genesymbol"="human_symbol")) %>%
        mutate(source_genesymbol = murine_symbol, .keep = "unused") %>%
        filter(!is.na(target_genesymbol) | !is.na(source_genesymbol)) %>%
        filter(!is.na(target_genesymbol)) %>%
        filter(!is.na(source_genesymbol))
}



#' Function to Tranform LIANA aggregated results
#' @param liana_agg a long tibble with liana_aggregate results
#'
#' @details converts `liana_aggregate` output to long with rank as predictor:
#'    source              target                   method_name     predictor
#     <chr>                <chr>                     <chr>            <dbl>
# Presomitic.mesoderm Anterior.somitic.tissues aggregate_rank          1
liana_agg_to_long <- function(liana_agg){
    liana_agg %>%
        mutate(across(c(source, target), ~str_replace_all(.x, " ", "."))) %>%
        dplyr::select(source, target, ends_with("rank"), -c(mean_rank, median_rank)) %>%
        mutate(aggregate_rank = min_rank(aggregate_rank)) %>% # convert to an integer rank
        pivot_longer(-c(source,target),
                     names_to = "method_name",
                     values_to = "predictor")
}




#' @title Calculate AUROC and PRROC from rank-adt tibbles- `prepare_for_roc`-formated
#' `adt_rank` elements
#'
#' @param df run_benchmark roc column provided as input
#' @param downsampling logical flag indicating if the number of Negatives
#'    should be downsampled to the number of Positives
#' @param times integer showing the number of downsampling
#' @param curve whether to return a Precision-Recall Curve ("PR") or ROC ("ROC")
#' @param seed An integer to set the RNG state for random number generation. Use
#'    NULL for random number generation.
#' @param source_name feature name (e.g. interaction, TF, etc)
#' @param auc_only whether to return auc only (excl th for sub_res/subsampling results)
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
                      seed = 1234,
                      source_name,
                      auc_only = FALSE,
                      ...){
    set.seed(seed)

    if(curve=="PR"){
        res_col_1 <- "precision"
        res_col_2 <- "recall"
        curve_fun = yardstick::pr_curve
        auc_fun = yardstick::pr_auc
    } else if(curve=="ROC"){
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
    feature_coverage = length(unique(df %>% pluck(source_name)))

    if(downsampling == TRUE){
        num_tp = nrow(cp)

        res = map_df(seq(from=1, to=times, by=1), function(i) {
            df_sub = sample_n(cn, num_tp, replace=TRUE) %>%
                bind_rows(cp)

            r_sub = df_sub %>%
                curve_fun(.data$response, .data$predictor, ...)

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
            curve_fun(.data$response, .data$predictor, ...)
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

    if(auc_only){
        res %<>%
            select(-c(th, !!res_col_1, !!res_col_2)) %>%
            distinct()
    }

    return(res)
}


#' Run Fischer's exact test on a contigency_table
#' @param cont_tab
enrich3 <- function(cont_tab){
    cont_table <- cont_tab %>%
        as.data.frame() %>%
        column_to_rownames("localisation") %>%
        t()

    result <- fisher.test(cont_table)
    tibble(pval = last(result$p.value), odds_ratio = result$estimate)
}



#' @title Recode dataset names
#' @param dataset - vector /w dataset names
recode_datasets <- function(datasets){
    dplyr::recode(datasets,
                  # CITE-Seq
                  "10k_malt" = "10kMALT",
                  "10k_pbmcs" = "10kPBMCs",
                  "5k_pbmcs" = "5kPBMCs ",
                  "5k_pbmcs_nextgem" = "5kPBMCs (nextgem)",
                  "cmbcs" = "3kCBMCs",
                  "spleen_lymph_101" = "SLN111",
                  "spleen_lymph_206" = "SLN208",

                  # Mouse Brain Cortex visium
                  "anterior1" = "Cortex Anterior 1",
                  "anterior2" = "Cortex Anterior 2",
                  "posterior1" = "Cortex Posterior 1",
                  "posterior2" = "Cortex Posterior 2",

                  "1142243F" = "TNBC1 (1142243F)",
                  "1160920F" = "TNBC2 (1160920F)",
                  "CID4290" = "ER1 (CID4290)",
                  "CID4465" = "TNBC3 (CID4465)",
                  "CID4535" = "ER2 (CID4535)",
                  "CID44971" = "TNBC4 (CID44971)"
                  )
}


#' Helper function to produce AUROC heatmap
#' @param roc_tibble Tibble with calculated AUROC/PRROC
#' @param curve type of curve `roc` or `prc`
#' @param mat_only whether to return only the auc_mat used to build the heatmap
#' in long format
#' @inheritDotParams ComplexHeatmap::Heatmap
#'
#' @return returns an AUROC or Precision-Recall AUC heatmap
#' @import ComplexHeatmap ggplot2 viridis
get_auroc_heat <- function(roc_tibble,
                           curve,
                           mat_only = FALSE,
                           auc_min = NULL,
                           auc_max = NULL,
                           ...){

    auc_df <- roc_tibble %>%
        dplyr::select(dataset, method_name, !!curve) %>%
        unnest(!!curve)

    auc_min %<>% `%||%` (min(auc_df$auc))
    auc_max %<>% `%||%` (max(auc_df$auc))

    auc_mat <- auc_df %>%
        dplyr::select(dataset, method_name, auc) %>%
        distinct() %>%
        mutate(method_name = gsub("\\..*", "", method_name)) %>%
        mutate(method_name = recode_methods(method_name)) %>%
        mutate(dataset = recode_datasets(dataset)) %>%
        pivot_wider(names_from = method_name, values_from = auc) %>%
        as.data.frame() %>%
        column_to_rownames("dataset") %>%
        as.matrix()

    if(mat_only){
        auc_mat %<>%
            as.data.frame() %>%
            pivot_longer(
                cols = everything(),
                names_to = "method",
                values_to = "estimate")

        return(auc_mat)
    }

    ComplexHeatmap::Heatmap(
        auc_mat,
        col = circlize::colorRamp2(c(auc_min,auc_max),
                                   viridis::cividis(2)),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        heatmap_legend_param =
            list(title = str_glue("AU{str_to_upper(curve)}",
                                  fontsize = 20,
                                  grid_height = unit(21, "mm"),
                                  grid_width = unit(21, "mm"))),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid::grid.text(sprintf("%.2f", auc_mat[i, j]),
                            x, y,
                            gp = grid::gpar(fontsize = 18,
                                            fontface = "bold",
                                            col = "white"))
        },
        ...
    )
}



#' Helper function to extend liana aggragate to e.g. filters of certain methods
#' @param liana_res in list form (not aggragated)
#' @param filt_de_pvals whether to filter differential expression p-values (e.g. Connectome)
#' @param de_thresh differential gene expression threshold
#' @param filt_outs whether to filter end output (CellChat, CellPhoneDB, SCA)
#' @param pval_thresh permutation methods output p-value threshold
#' @param sca_thresh SingleCellSignalR threshold
#' @inheritDotParams liana_aggregate
liana_aggregate_enh <- function(liana_res,
                                filt_de_pvals = TRUE,
                                de_thresh = 0.05,
                                filt_outs = FALSE,
                                pval_thresh = 0.05,
                                sca_thresh = 0.5,
                                ...){
    if(filt_de_pvals){
        # Filter Connectome genes below 0.05
        liana_res$call_connectome %<>%
            filter(p_val_adj.lig <= de_thresh) %>%
            filter(p_val_adj.rec <= de_thresh)

    }

    # Filter according to end threshold
    if(filt_outs){
        ## CellChat
        liana_res$cellchat %<>% filter(pval <= pval_thresh)
        ## CellPhoneDB/Squidpy
        liana_res$squidpy %<>% filter(pvalue <= pval_thresh)
        ## sca
        liana_res$call_sca %<>% filter(LRscore >= sca_thresh)
    }

    # aggregate liana results
    liana_res %<>% liana_aggregate(...,
                                   cap = 500000)

    return(liana_res)
}
