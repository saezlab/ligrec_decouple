#### Run LIANA with all original methods (with defaults args) + multiple resources
require(tidyverse)
require(magrittr)
require(Seurat)
require(liana)

source("src/eval_utils.R")

# ER
brca_subtype <- "ER"
liana_res <- readRDS("data/output/comparison_out/BRCA_ER_liana_res.RDS")
lr_er <- liana_res %>%
    transpose() %>%
    pluck("OmniPath") %>%
    liana_aggregate_enh(filt_de_pvals = TRUE,
                        de_thresh = 0.05, # we only filter connectome DEs
                        filt_outs = FALSE,
                        pval_thresh = 1,
                        sca_thresh = 0,
                        .score_mode = liana:::.score_specs)


# TNBC
brca_subtype <- "TNBC"
liana_res <- readRDS("data/output/comparison_out/BRCA_TNBC_liana_res.RDS")
lr_tbnc <- liana_res %>%
    transpose() %>%
    pluck("OmniPath") %>%
    liana_aggregate_enh(filt_de_pvals = TRUE,
                        de_thresh = 0.05, # we only filter connectome DEs
                        filt_outs = FALSE,
                        pval_thresh = 1,
                        sca_thresh = 0,
                        .score_mode = liana:::.score_specs)



# HER2
brca_subtype <- "HER2"
liana_res <- readRDS("data/output/comparison_out/BRCA_HER2_liana_res.RDS")
lr_her2 <- liana_res %>%
    transpose() %>%
    pluck("OmniPath") %>%
    liana_aggregate_enh(filt_de_pvals = TRUE,
                        de_thresh = 0.05, # we only filter connectome DEs
                        filt_outs = FALSE,
                        pval_thresh = 1,
                        sca_thresh = 0,
                        .score_mode = liana:::.score_specs)
