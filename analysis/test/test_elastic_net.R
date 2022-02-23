require(tidyverse)
require(magrittr)
require(decoupleR)
require(Seurat)
require(SingleCellExperiment)

# load old results
tnbc_brca <- readRDS("data/output/cytosig_out/BRCA_TNBC_cytosig.RDS")


source("analysis/cytosig/cytosig_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

pseudobulk <- get_pseudobulk(readRDS("data/input/spatial/Wu_etal_2021_BRCA/deconv/TNBC_celltype_minor/TNBC_celltype_minor_seurat.RDS"),
                             expr_prop = 0.1,
                             assay = "RNA",
                             sum_count_thresh = 5)
saveRDS(pseudobulk, "~/Downloads/pseudobulk_test.RDS")

cytosig_net <- load_cytosig() %>%
    mutate(mor = mor*weight)



# Pseudobulk cytosig
pseudo_cytosig <- pseudobulk %>%
    mutate(cytosig_res = logcounts %>%
               map(function(logc){
                   run_mlm(
                       logc,
                       cytosig_net,
                       .source = "cytokine",
                       .target = "target",
                       .mor = "mor",
                       sparse = FALSE
                   ) %>%
                       # rename
                       select(cytokine=source,
                              NES=score,
                              p_value)
               })
    ) %>%
    select(celltype, cytosig_res)


### ELASTIC NET
# Pseudobulk cytosig
elastic_cytosig <- pseudobulk %>%
    mutate(cytosig_res = logcounts %>%
               map(function(logc){
                   run_mlm(
                       logc,
                       cytosig_net,
                       .source = "cytokine",
                       .target = "target",
                       .mor = "mor",
                       sparse = FALSE,
                       regularization = 'ridge'
                   ) %>%
                       # rename
                       select(cytokine=source,
                              NES=score,
                              p_value)
               })
    ) %>%
    select(celltype, cytosig_res)


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



cytosig_res <- elastic_cytosig %>%
    unnest(cytosig_res) %>%
    # keep only the most positive cytokines
    group_by(celltype) %>%
    top_n(n=5, wt=NES) %>%
    # join aliases
    left_join(alias_tib, by = "cytokine") %>%
    mutate(aliases = if_else(is.na(aliases),
                             cytokine,
                             aliases)) %>%
    select(-cytokine, cytokine=aliases) %>%
    select(cytokine, celltype, NES) %>%
    ungroup(celltype)

# load liana
liana_path <- readRDS(file.path("data/output/aggregates", str_glue("TNBC_independent_mixed_liana_res.RDS")))





