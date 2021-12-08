# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(patchwork)

source("analysis/cytosig/cytosig_src.R")
source("src/plot_utils.R")
source("src/eval_utils.R")
source("analysis/spatial/spatial_src.R")

# Function to Combine the plots
cyto_space_patch <- function(cytosig_p,
                             space_p,
                             path
                             ){
    cairo_pdf(path,
              height = 24,
              width = 22,
              family = 'DINPro')
    print((cytosig_p / space_p) +
              plot_layout(guides = 'keep', heights = c(1, 1)) +
              plot_annotation(tag_levels = 'A',
                              tag_suffix = ')') &
              theme(plot.tag = element_text(face = 'bold',
                                            size = 40)))
    dev.off()
}



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

#
# # Specificities Method specifics
# space_specs <- readRDS("data/output/spatial_out/all_fets_specs.RDS") %>%
#     get_spatial_boxplot()
# cytosig_specs <- plot_cytosig_aucs(.eval = "independent", score_mode = "specs")
#
# path <- file.path("figures",
#                   "SuppFigure20_Evals_Specific.RDS")
# cyto_space_patch(space_specs,
#                  cytosig_specs,
#                  path)
#
#
# Harmonized Plot
cytosig_harmonize_data <- get_cytosig_fets(
    inputpath = "data/output/eval_harmonize/cytosig_res_independent_comp.RDS"
    )
cytosig_harmonize <- cytosig_harmonize_data %>%
    get_eval_boxplot(eval_type = "cytosig")
space_harmonize <- readRDS("data/output/eval_harmonize/harmonize_lr_coloc.RDS") %>%
    get_eval_boxplot(eval_type = "space")

path <- file.path("figures",
                  "SuppFigure21_Harmonize_CytoSpace.RDS")
cyto_space_patch(cytosig_harmonize,
                 space_harmonize,
                 path)
