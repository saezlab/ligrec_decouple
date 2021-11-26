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
    print(cytosig_plot(.eval = .eval,
                 score_mode = score_mode))
    dev.off()
}

# Generate temp plots
comb_tibble %>%
    pmap(~print_cyto_plot(.x, .y))


liana_res <- readRDS("data/output/comparison_out/BRCA_ER_liana_res.RDS") %>%
    transpose() %>%
    pluck("OmniPath")
gc()

liana_agg <- liana_res %>%
    liana_aggregate_enh(filt_de_pvals = TRUE,
                        de_thresh = 0.05,
                        filt_outs = TRUE,
                        pval_thresh = 1,
                        sca_thresh = 0,
                        .eval = "independent")

agg_sample <- liana_agg %>%
    head(250000) %>%
    rowwise() %>%
    mutate(
        across(
            ends_with("rank"),
            # max is imputed by liana_aggregate by default
            # thus, we simply set it as NA if we don't want to consider it
            ~na_if(.x, cap))
        ) %>%
    mutate(
        across(
            ends_with("rank"),
            # max is imputed by liana_aggregate by default
            # thus, we simply set it as NA if we don't want to consider it
            ~na_if(.x, nrow(liana_res)))
    )

rowwise() %>%
    mutate(
        across(
            ends_with("rank"),
            # max is imputed by liana_aggregate by default
            # thus, we simply set it as NA if we don't want to consider it
            ~ifelse(.x==(nrow(liana_agg) || .x==cap), NA, .x)
        )) %>%
    ungroup()

#^ Try with max and compare
maxsp <- readRDS("data/output/brca_extracts/ER_max_specs_liana_res.RDS")

all_equal(liana_res, maxsp)

ind

indp <- readRDS("data/output/brca_extracts/ER_independent_specs_liana_res.RDS")


