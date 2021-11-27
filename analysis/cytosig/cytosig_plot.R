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


