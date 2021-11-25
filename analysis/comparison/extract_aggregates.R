# Script Used to Extract and Aggregate OmniPath results for each dataset
# used in the evaluations
library(tidyverse)
library(liana)
library(Seurat)
library(magrittr)

source("src/eval_utils.R")
source("src/plot_utils.R")

# Get Args from std in
args <- commandArgs(trailingOnly=TRUE)
dataset <- args[[2]] # also job name # (e.g. "HER2")
liana_path <- args[[3]] # path to liana (e.g. "data/output/comparison_out/BRCA_HER2_liana_res.RDS")
eval <- args[[4]] # eval type ("intersect", "independent", "max")

# Read Liana results
liana_res <- readRDS(liana_path) %>%
    # Only keep OP
    transpose() %>%
    pluck("OmniPath")
gc()

# I) Specifics
liana_specs_agg <- liana_res %>%
    liana_aggregate_enh(
        pval_thresh = 1,
        sca_thresh = 0,
        de_thresh = 0.05,
        .score_mode = liana:::.score_specs,
        .eval = eval
        )
saveRDS(liana_specs_agg, file.path("data/output", str_glue("{dataset}_{eval}_specs_liana_res.RDS")))
rm(liana_specs_agg)
gc()

# II) Housekeeping
liana_house_agg <- liana_res %>%
    liana_aggregate_enh(
        pval_thresh = 0.05,
        sca_thresh = 0,
        de_thresh = 0.05,
        .score_mode = liana:::.score_housekeep,
        .eval = eval
)
saveRDS(liana_house_agg, file.path("data/output", str_glue("{dataset}_{eval}_specs_liana_res.RDS")))
rm(liana_house_agg)
gc()


# III) Comp/Mixed Scores
liana_mixed_agg <- liana_res %>%
    liana_aggregate_enh(
        pval_thresh = 0.05,
        sca_thresh = 0,
        de_thresh = 0.05,
        .score_mode = .score_comp,
        .eval = eval
)
saveRDS(liana_mixed_agg, file.path("data/output", str_glue("{dataset}_{eval}_specs_liana_res.RDS")))
rm(liana_mixed_agg)
gc()



