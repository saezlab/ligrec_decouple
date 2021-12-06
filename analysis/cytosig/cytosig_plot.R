# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)
require(SingleCellExperiment)
require(decoupleR)

source("analysis/cytosig/cytosig_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

# Loop over all plots
eval_vec <- c("independent", "max", "intersect")
score_mode_vec <- c("mixed", "specs", "house")

comb_tibble <- expand_grid(eval_vec, score_mode_vec)

print_cyto_plot <- function(.eval,
                            score_mode){
    cairo_pdf(file.path("data/output/temp/",
                        str_glue("cyto_{.eval}_{score_mode}_.pdf")),
              height = 12,
              width = 16,
              family = 'DINPro')
    print(plot_cytosig_aucs(.eval = .eval,
                 score_mode = score_mode))
    dev.off()
}


#### Purgatory ----
enrich3 <- function(cont_tab){
    cont_table <- cont_tab %>%
        as.data.frame() %>%
        column_to_rownames("response") %>%
        t()

    result <- fisher.test(cont_table)
    tibble(pval = last(result$p.value), odds_ratio = result$estimate)
}



xx <- readRDS("data/output/cytosig_out/cytosig_res_independent_mixed.RDS")
cytolr <- xx$cytosig_res[[1]] %>% unnest(cyto_liana) %>%
    select(-c(roc, prc, corr))

n_rank = 100

ranks_counted <- cytolr %>%
    group_by(method_name) %>%
    mutate(predictor = min_rank(predictor*-1)) %>%
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
                       cont_tab %>% enrich3
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
    arrange(desc(enrichment))
fet_results


