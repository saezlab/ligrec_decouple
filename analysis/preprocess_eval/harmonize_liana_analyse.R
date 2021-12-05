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

# 1) Cytosig ----
.eval = "independent"
# Customize tibble to harmonized results
path_tibble <-
    tibble(dataset = c("ER",
                       "HER2",
                       "TNBC"),
           # BRCA seurat objects used throughout the manuscript
           seurat_path = c("data/input/spatial/Wu_etal_2021_BRCA/deconv/ER_celltype_minor/ER_celltype_minor_seurat.RDS",
                           "data/input/spatial/Wu_etal_2021_BRCA/deconv/HER2_celltype_minor/HER2_celltype_minor_seurat.RDS",
                           "data/input/spatial/Wu_etal_2021_BRCA/deconv/TNBC_celltype_minor/TNBC_celltype_minor_seurat.RDS"
           ),
           # liana results generated from extract_evals.sh
           liana_path = c(file.path("data/output/eval_harmonize", str_glue("BRCA_ER_liana_agg.RDS")),
                          file.path("data/output/eval_harmonize", str_glue("BRCA_HER2_liana_agg.RDS")),
                          file.path("data/output/eval_harmonize", str_glue("BRCA_TNBC_liana_agg.RDS"))
           ))

# Run Cytosig pipe
cytosig_eval_wrap(.eval = .eval,
                  score_mode = "comp",
                  generate = FALSE,
                  path_tibble = path_tibble,
                  outpath = "data/output/eval_harmonize")

# Plot results
plot_cytosig_aucs(inputpath = "data/output/eval_harmonize/cytosig_res_independent_comp.RDS")



# 2) Spatial ----
corr_thresh <- 1.645
n_ranks = c(50, 100,
            500, 1000,
            2500, 5000,
            10000)

visium_dict <- list("ER",
                    "TNBC",
                    "brain"
                   )

# Tibble with Paths
coloc_tibble <- tibble::tribble(~condition, ~liana_agg_path, ~spatial_corr_path,
                                "ER+ BRCA",
                                "data/output/eval_harmonize/BRCA_ER_liana_agg.RDS",
                                "data/output/spatial_out/deconv_summ/ER_coloc_corr.RDS",
                                "TNBC",
                                "data/output/eval_harmonize/BRCA_TNBC_liana_agg.RDS",
                                "data/output/spatial_out/deconv_summ/TNBC_coloc_corr.RDS",
                                "Brain Cortex",
                                "data/output/eval_harmonize/brain_liana_agg.RDS",
                                "data/output/spatial_out/deconv_summ/brain_coloc_corr.RDS"
                                )

# Tibble with FET results
coloc_tibble %<>%
    pmap(function(condition, liana_agg_path, spatial_corr_path){
        get_lr_colocalized(liana_agg_path,
                           spatial_corr_path,
                           condition = condition,
                           corr_thresh = 1.645,
                           n_ranks = c(50, 100, 500, 1000,
                                       2500, 5000, 10000))
    }) %>%
    bind_rows()
saveRDS(coloc_tibble,"data/output/eval_harmonize/harmonize_lr_coloc.RDS")

# Spatial Plot
coloc_tibble %>%
    get_spatial_boxplot()
