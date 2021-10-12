#' Function to reshape estimate to long format
#' @param estimate_mat diagonal df/matrix with colocalisation estimates
#' @param z_scale whether to z-tranform the estimate column
#'
#' @details converts the diagonal matrix/df to a long format tibble with
#' celltype1, celltype2, and estimate (e.g. correlation/z-scores) columns
reshape_coloc_estimate <- function(estimate_mat, z_scale = FALSE){
    estimate_mat %>%
        reshape2::melt() %>%
        as_tibble() %>%
        dplyr::rename(celltype1 = Var1,
                      celltype2 = Var2,
                      estimate = value) %>%
        # FILTER AUTOCRINE
        filter(celltype1!=celltype2) %>%
        mutate(estimate = ifelse(rep(z_scale, length(estimate)),
                                 scale(estimate),
                                 estimate))
}



#' Function to perform FET on  colocalized vs not_colocalized interactions
#'
#' @param coloc_estimate tibble with colocalization metric for each pair of cell types
#' in long format: cols = c(celltype1, celltype2, estimate)
#' @param liana_format liana in long format with cols:
#' source target method_name predictor(rank)
#' @param n_rank number of ranks to be considered for the top interactions
#' @param arb_thresh a threshold for the co-localization estimate
#'  (i.e. correlation/z-score above or below X)
#'
#' @details What I do is that, I compare (FET) the proportion of colocalized cell types
#' in the top X interactions to the number of colocalized cell types in all interactions.
#' (top vs total)
#' With the assumption that a good method should have a higher proportion of
#' colocalized cell types in its top rank hits.
#'
#' Cont table example:
#'>    localisation       total   top
#'>       <chr>          <int>    <int>
#'>  1 not_colocalized   299868    748
#'>  2 colocalized        20834    159
run_coloc_fet <- function(coloc_estimate, liana_format, n_rank, arb_thresh){

    # Assign colocalisation to liana results in long according to a threshold
    liana_loc <- liana_format %>%
        left_join(coloc_estimate, by = c("source"="celltype1",
                                          "target"="celltype2"))  %>%
        # FILTER AUTOCRINE
        filter(source!=target) %>%
        dplyr::mutate(localisation = case_when(estimate >= arb_thresh ~ "colocalized",
                                               estimate < arb_thresh ~ "not_colocalized"
        )) %>% ungroup()

    # count total vs top colocalized
    ranks_counted <- liana_loc %>%
        # count total (Null)
        group_by(method_name, localisation) %>%
        mutate(total = n())  %>%
        # count in x rank (Alt)
        filter(predictor <= n_rank) %>%
        group_by(method_name, localisation) %>%
        mutate(top = n()) %>%
        dplyr::select(method_name, localisation, total, top) %>%
        distinct()

    # contingency tables for colocalized vs not by method
    cont_tabs <- ranks_counted %>%
        group_by(method_name) %>%
        group_nest(.key = "contigency_tables")
    cont_tabs$contigency_tables[[1]]

    # Fischer's exact test on co-localized in top vs total interactions
    fet_results <- cont_tabs %>%
        # FET
        mutate(fet_res = contigency_tables %>%
                   map(function(cont_tab) cont_tab %>% enrich3)) %>%
        # unnest fet results
        dplyr::select(-contigency_tables) %>%
        unnest(fet_res) %>%
        # modify results
        mutate(
            padj = p.adjust(pval, method = "fdr"),
            enrichment = ifelse(
                odds_ratio < 1,
                -1 / odds_ratio,
                odds_ratio
            ) %>% unname
        ) %>%
        arrange(desc(enrichment))
    fet_results

}
