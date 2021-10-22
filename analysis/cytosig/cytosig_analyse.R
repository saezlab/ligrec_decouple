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

# get cytosig
cytosig_net <- load_cytosig()

# Load Seurat Object
seurat_object <- readRDS("data/input/spatial/Wu_etal_2021_BRCA/deconv/ER_celltype_minor/ER_celltype_minor_seurat.RDS")

# read liana
liana_res <- readRDS("data/input/spatial/Wu_etal_2021_BRCA/deconv/ER_celltype_minor/ER_celltype_minor_liana_res.RDS")

cyto_res <- run_cytosig_eval(seurat_object, liana_res,
                             cytosig_net, z_scale = TRUE,
                             expr_prop = 0.1, assay = "RNA",
                             sum_count_thresh = 5, NES_thresh = 1.645)




# scale
z_scale = TRUE

# Load Seurat Object
seurat_object <- readRDS("data/input/spatial/Wu_etal_2021_BRCA/deconv/ER_celltype_minor/ER_celltype_minor_seurat.RDS")

# read liana
liana_res <- readRDS("data/input/spatial/Wu_etal_2021_BRCA/deconv/ER_celltype_minor/ER_celltype_minor_liana_res.RDS") %>%
    liana_aggregate()


# get cytosig
cytosig_net <- load_cytosig()

# get pseudobulk
pseudo <- get_pseudobulk(seurat_object,
                         assay = "RNA",
                         expr_prop = 0.1,
                         sum_count_thresh = 5)

# Cytokine Activity Enrichment
pseudo_cytosig <- pseudo %>%
    mutate(cytosig_res = logcounts %>%
               map(function(logc){
                   run_wmean(
                       logc,
                       cytosig_net,
                       .source = "cytokine",
                       .target = "target",
                       .mor = "mor",
                       .likelihood = "weight",
                       times = 1000,
                       seed = 1234,
                       sparse = TRUE,
                       randomize_type = "cols_independently") %>%
                       # keep only norm weighted mean
                       filter(statistic == "norm_wmean") %>%
                       # rename
                       select(cytokine=source,
                              NES=score,
                              p_value) %>%
                       # correct p
                       mutate(adj_pvalue = p.adjust(p_value))
               }))
saveRDS(pseudo_cytosig, "data/output/cytosig_out/ER_BRCA_cytosig.RDS")

### Join LIANA to cytosig
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

# cytosig results format
cytosig_res <- pseudo_cytosig %>%
    select(celltype, cytosig_res) %>%
    unnest(cytosig_res) %>%
    left_join(alias_tib, by = "cytokine") %>%
    mutate(aliases = if_else(is.na(aliases),
                             cytokine,
                             aliases)) %>%
    select(-cytokine, cytokine=aliases) %>%
    select(cytokine, celltype, NES) %>% # exludes p-vals
    { if(z_scale) z_transform_nes(.)  else . } # cluster-specific or not


z_transform_nes <- function(df){
    df %>%
        group_by(cytokine) %>%
        column_to_rownames("cytokine") %>%
        mutate(NES = scale(NES)) %>%
        unnest(NES)
}


# Format and join cytokine activities and prep for AUC
liana_cytosig <- liana_res %>%
    # keep only interactions in which the ligand is a cytokine present in cytosig
    filter(ligand %in% unique(cytosig_res$cytokine)) %>%
    dplyr::select(ligand, target, ends_with("rank"), -c(mean_rank, median_rank)) %>%
    mutate(aggregate_rank = min_rank(aggregate_rank)) %>%
    # here we join by ligand from liana and cytokine in cytosig,
    # as well as target cluster from liana and the celltype for which we predicted
    # the ligand activities
    left_join(cytosig_res, by=c("ligand"="cytokine", "target"="celltype")) %>%
    pivot_longer(ends_with("rank"), names_to = "method_name", values_to = "predictor") %>%
    mutate(predictor = predictor*-1) %>%
    unite(ligand, target, col = "cytokine_in_target")


# Prep for ROC
cytosig_eval <- liana_cytosig %>%
    mutate(response = if_else(NES > 1.7,
                              1,
                              0)) %>%
    mutate(response = factor(response, levels = c(1, 0))) # first level is the truth

# ROC and Correlations
cytosig_eval %<>%
    group_by(method_name) %>%
    group_nest(.key = "cyto_liana") %>%
    mutate(roc = map(cyto_liana,
                     function(df) calc_curve(df,
                                             curve="ROC",
                                             downsampling = FALSE,
                                             times = 100,
                                             source_name = "cytokine_in_target",
                                             auc_only = FALSE
                     )
    )) %>%
    mutate(prc = map(cyto_liana,
                     function(df) calc_curve(df,
                                             curve="PR",
                                             downsampling = TRUE,
                                             times = 1000,
                                             source_name = "cytokine_in_target",
                                             auc_only = FALSE))) %>%
    mutate(corr = cyto_liana %>%
               map(function(df){
                   cor.test(df[["NES"]],
                            df[["predictor"]],
                            method="spearman",
                            exact = FALSE
                   ) %>% broom::tidy()
               }))
# ^ This into pipe


