# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(patchwork)

source("analysis/cytosig/cytosig_src.R")
source("src/plot_utils.R")
source("src/eval_utils.R")
source("analysis/spatial/spatial_src.R")

# Mixed/Comp Method specifics
cytosig_mixed_data <- get_cytosig_fets(.eval = "independent",
                                       score_mode = "mixed")
cytosig_mixed <- cytosig_mixed_data %>%
    get_eval_boxplot(eval_type = "cytosig")
#plot
space_mixed <- readRDS("data/output/spatial_out/all_fets_mixed.RDS") %>%
    get_eval_boxplot(eval_type="space")
path <- file.path("figures",
                  "Figure5_Evals_Composite.RDS")
cyto_space_patch(cytosig_mixed,
                 space_mixed,
                 path)


