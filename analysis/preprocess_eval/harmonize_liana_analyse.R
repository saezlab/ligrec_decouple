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
source("src/eval_utils.R")
source("src/plot_utils.R")

### Cytosig ----
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

### Spatial ----

