# Script to Re-do the Cytosig and Visium-Spatial evaluations
# with harmonized filtering

# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)
require(SingleCellExperiment)
require(decoupleR)

source("analysis/cytosig/cytosig_src.R")
source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

# Cytosig: add aliases (as in OP) and cytokine family members if appropriate
aliases_and_families <- list("CD40L" = "CD40LG",
                             "GSFC" = "CSF3",
                             "IFN1" = c("IFNA1", "IFNA2", "IFNA10",
                                        "IFNA7", "IFNA21", "IFNA5",
                                        "IFNA14", "IFNA17", "IFNA6",
                                        "IFNA4",  "IFNA16", "IFNA8" ),
                             "IFNL" = c("IFNL1", "IFNL2", "IFNL3", "IFNL4"),
                             "IL12" = c("IL12A", "IL12B"),
                             "IL36" = c("IL36A", "IL36B", "IL36G", "IL36RN"),
                             "MCSF" = "CSF1",
                             "TNFA" = c("TNF", "TNFA"),
                             "TRAIL" = c("TNFSF10"),
                             "TWEAK" = "TNFSF12")
alias_tib <- tibble(cytokine = names(aliases_and_families),
                    aliases = aliases_and_families) %>%
    unnest(aliases)

#
z_scale = FALSE
cytosig_res <- readRDS("data/output/cytosig_out/BRCA_HER2_cytosig.RDS") %>%
    unnest(cytosig_res) %>%
    left_join(alias_tib, by = "cytokine") %>%
    mutate(aliases = if_else(is.na(aliases),
                             cytokine,
                             aliases)) %>%
    select(-cytokine, cytokine=aliases) %>%
    select(cytokine, celltype, NES, p_value) %>%
    { if(z_scale) z_transform_nes(.)  else . }

#
cytotalk_test <-# readRDS("data/output/cytosig_out/HER") %>%
    #mutate(crosstalk_rank = min_rank(-crosstalk_score))
    readRDS("data/output/aggregates/HER2_independent_mixed_liana_res.RDS")


# Format and join cytokine activities and prep for AUC
liana_cytosig <- cytotalk_test %>%
    # keep only interactions in which the ligand is a cytokine present in cytosig
    filter(ligand %in% unique(cytosig_res$cytokine)) %>%
    dplyr::select(ligand, target, ends_with("rank")) %>%
    # here we join by ligand from liana and cytokine in cytosig,
    # as well as target cluster from liana and the celltype for which we predicted
    # the ligand activities
    left_join(cytosig_res, by=c("ligand"="cytokine", "target"="celltype")) %>%
    pivot_longer(ends_with("rank"), names_to = "method_name", values_to = "predictor") %>%
    mutate(predictor = predictor*-1) %>% # convert rankings
    unite(ligand, target, col = "cytokine_in_target")


# Prep for ROC
NES_thresh = 1.645
cytosig_eval <- liana_cytosig %>%
    # correct p
    mutate(adj_pvalue = p.adjust(p_value, "fdr")) %>%
    mutate(response = if_else(NES > NES_thresh & adj_pvalue <=0.05,
                              1,
                              0)) %>%
    mutate(response = factor(response, levels = c(1, 0))) # first level is the truth
print(cytosig_eval %>% arrange(NES))



#
cytosig_eval


message("Calculating AUCs")
# ROC and Correlations
cytosig_eval %<>%
    group_by(method_name) %>%
    group_nest(.key = "cyto_liana") %>%
    mutate(roc = map(cyto_liana,
                     function(df) calc_curve(df,
                                             curve="ROC",
                                             downsampling = FALSE,
                                             source_name = "cytokine_in_target",
                                             auc_only = TRUE)
    )) %>%
    mutate(prc = map(cyto_liana,
                     function(df) calc_curve(df,
                                             curve="PR",
                                             downsampling = TRUE,
                                             times = 100,
                                             source_name = "cytokine_in_target",
                                             auc_only = TRUE))) %>%
    mutate(corr = cyto_liana %>%
               map(function(df){
                   cor.test(df[["NES"]],
                            df[["predictor"]],
                            method="spearman",
                            exact = FALSE) %>%
                       broom::tidy()
               }))
cytosig_eval$roc[[1]]


cytolr <- cytosig_eval %>%
    select(-c(roc, prc, corr))  %>%
    unnest(cyto_liana) %>%
    group_by(method_name) %>%
    mutate(rank = predictor*-1) %>%
    mutate(predictor = min_rank(predictor*-1))

max_ranks <- cytolr %>%
    select(method_name, predictor) %>%
    group_by(method_name) %>%
    na.omit() %>%
    summarise(max_rank = max(predictor))

# FET ----
n_ranks = c(100, 250, 500, 1000,
            2500, 5000, 10000)
map(n_ranks, function(n_rank){
    # Count total vs top
    ranks_counted <- cytolr %>%
        group_by(method_name, response) %>%
        mutate(total = n())  %>%
        # count in x rank (Alt)
        filter(predictor <= n_rank) %>%
        group_by(method_name, response) %>%
        mutate(top = n()) %>%
        dplyr::select(method_name, response, total, top) %>%
        distinct() %>%
        ungroup()

    # contingency tables for colocalized vs not by method
    cont_tabs <- ranks_counted %>%
        group_by(method_name) %>%
        group_nest(.key = "contigency_tables")

    # Fischer's exact test on co-localized in top vs total interactions
    fet_results <- cont_tabs %>%
        # FET
        mutate(fet_res = contigency_tables %>%
                   map(function(cont_tab){
                       # handle case where no colocolalized are in top X ranks
                       if(nrow(cont_tab) == 1){
                           if(cont_tab %>% pluck("response")==0){
                               message("Only negative class is present")
                               tibble(pval=1,
                                      odds_ratio=-9999)
                           } else{
                               stop("Only positive class is present!!!")
                           }
                       } else{ # run enrichment
                           cont_tab %>%
                               arrange(desc(response)) %>%
                               enrich3(., "response")
                       }
                   })) %>%
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
        arrange(desc(enrichment)) %>%
        mutate(n_rank = n_rank)
}) %>%
    bind_rows() -> xxx

max_ranks
xxx

xxx2 <- xxx %>%
    left_join(max_ranks) %>%
    mutate(n_rank = ifelse(n_rank > max_rank,max_rank, n_rank)) %>%
    mutate(method_name = gsub("\\..*","", method_name)) %>%
    mutate(method_name = recode_methods(method_name))
