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


seurat_object <- readRDS("data/input/spatial/Wu_etal_2021_BRCA/deconv/HER2_celltype_minor/HER2_celltype_minor_seurat.RDS")
pseudo <- get_pseudobulk(seurat_object,
                         assay = "RNA",
                         expr_prop = 0.1,
                         sum_count_thresh = 5)
cytosig_net <- load_cytosig(n_gene = 500)

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
                       randomize_type = "cols_independently"
                   ) %>%
                       # keep only norm weighted mean
                       filter(statistic == "norm_wmean") %>%
                       # rename
                       select(cytokine=source,
                              NES=score,
                              p_value)
               }))

## WMEAN
pseudo_cytosig %>%
    select(celltype, cytosig_res) %>%
    unnest(cytosig_res) %>%
    mutate(p_adj = p.adjust(p_value, "fdr")) %>%
    arrange(p_adj) %>%
    mutate(flag = p_adj <= 0.05 & NES > 0) %>%
    group_by(flag) %>%
    summarise(n())


### test MLM ----
seurat_object <- readRDS("data/input/spatial/Wu_etal_2021_BRCA/deconv/HER2_celltype_minor/HER2_celltype_minor_seurat.RDS")
pseudo <- get_pseudobulk(seurat_object,
                         assay = "RNA",
                         expr_prop = 0.1,
                         sum_count_thresh = 5)
cytosig_net <- load_cytosig(n_gene = 500)

pseudo_cytosig_mlm <- pseudo %>%
    mutate(cytosig_res = logcounts %>%
               map(function(logc){
                   run_mlm(
                       logc,
                       cytosig_net,
                       .source = "cytokine",
                       .target = "target",
                       .mor = "mor",
                       .likelihood = "weight",
                       sparse = FALSE
                   ) %>%
                       # rename
                       select(cytokine=source,
                              NES=score,
                              p_value)
               }))
pseudo_cytosig_mlm %>%
    select(celltype, cytosig_res) %>%
    unnest(cytosig_res) %>%
    mutate(p_adj = p.adjust(p_value, "fdr")) %>%
    arrange(p_adj) %>%
    mutate(flag = p_adj <= 0.05 & NES > 0) %>%
    group_by(flag) %>%
    summarise(n())


#### ----
# Another object
seurat_object <- readRDS("data/input/spatial/Wu_etal_2021_BRCA/deconv/ER_celltype_minor/ER_celltype_minor_seurat.RDS")
pseudo <- get_pseudobulk(seurat_object,
                         assay = "RNA",
                         expr_prop = 0.1,
                         sum_count_thresh = 5)
cytosig_net <- load_cytosig(n_gene = 500)


pseudo_cytosig <- pseudo %>%
    mutate(cytosig_res = logcounts %>%
               map(function(logc){
                   run_mlm(
                       logc,
                       cytosig_net,
                       .source = "cytokine",
                       .target = "target",
                       .mor = "mor",
                       .likelihood = "weight",
                       sparse = FALSE
                   ) %>%
                       # rename
                       select(cytokine=source,
                              NES=score,
                              p_value)
               }))
pseudo_cytosig %>%
    select(celltype, cytosig_res) %>%
    unnest(cytosig_res) %>%
    mutate(p_adj = p.adjust(p_value, "fdr")) %>%
    arrange(p_adj) %>%
    mutate(flag = p_adj <= 0.05 & NES > 0) %>%
    group_by(flag) %>%
    summarise(n())

pseudo_cytosig %<>%
    select(celltype, cytosig_res) %>%
    unnest(cytosig_res) %>%
    mutate(p_adj = p.adjust(p_value, "fdr")) %>%
    arrange(p_adj)   %>%
    mutate(response = if_else(NES > 0 & p_adj <=0.05,
                              1,
                              0))


