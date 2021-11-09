source("analysis/cytosig/cytosig_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

require(tidyverse)
require(liana)
library(decoupleR)
require(Seurat)
require(SingleCellExperiment)

inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")

mat <- file.path(inputs_dir, "input-expr_matrix.rds") %>%
    readRDS() %>%
    dplyr::glimpse()

cytosig_net <- load_cytosig(cytosig_path = "~/Repos/ligrec_decouple/data/input/cytosig/cytosig_signature_centroid.csv")

seurat_object <- readRDS("data/input/spatial/Wu_etal_2021_BRCA/deconv/HER2_celltype_minor/HER2_celltype_minor_seurat.RDS")
pseudo <- get_pseudobulk(seurat_object,
                         assay = "RNA",
                         expr_prop = 0.1,
                         sum_count_thresh = 10)

decres <- run_wmean(
    mat = pseudo$logcounts[[1]],
    network = cytosig_net,
    .source = "cytokine",
    .target = "target",
    .mor = "mor",
    .likelihood = "weight",
    # randomize_type = "cols_independently",
    times = 1000,
    seed = 1234,
    sparse = FALSE
    ) %>%
    filter(statistic == "norm_wmean")
