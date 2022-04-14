## Script to generate the comparison data presented in the Manuscript
library(tidyverse)
library(magrittr)
library(liana)
library(RColorBrewer)
library(pheatmap)
library(proxy)
library(UpSetR)
library(grid)
library(ComplexHeatmap)
library(patchwork)
library(ggplotify)

source("analysis/comparison/comparison_utils.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

#' @title Recode method names
#' @param dataset - vector /w dataset names
recode_datasets <- function(datasets){
    dplyr::recode(datasets,
                  !!!as.list(.dataset_keys)
    )
}
comparison_out <- "data/output/comparison_out/"


## I. Comp_n ----
comp_tibble <- comp_summ_plot(pattern = "comp_n",
                              comparison_out = comparison_out,
                              box_name = "Figure5.pdf",
                              heat_name = "SuppFig16_mixed_n_JI_heat.pdf")
comp_tibble %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))
# Note CBMC has ~identical Jaccard for different resources and methods

# Save Source Data
comp_tibble


## II. Comp_frac ----
comp_tibble_frac <- comp_summ_plot(pattern = "comp_frac",
                                   comparison_out = comparison_out,
                                   box_name = "SuppFig_12_Composite_frac.pdf",
                                   heat_name = "SuppFig17_mixed_frac_JI_heat.pdf")
comp_tibble_frac %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))


## III. Specs_n ----
n_tibble <- comp_summ_plot(pattern = "specs_n",
                           comparison_out = comparison_out,
                           box_name = "SuppFig_18_Specificity_n.pdf",
                           heat_name = "SuppFig19_specs_n_JI_heat.pdf")
n_tibble %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))


## IV. Specs_FRAC ----
# frac_tibble <- comp_summ_plot(pattern = "specs_frac",
#                               comparison_out = comparison_out,
#                               box_name = "SuppFig_12_Specificity_frac.pdf",
#                               heat_name = "SuppFig12_specs_frac_JI_heat.pdf")
# frac_tibble %>%
#     group_by(entity) %>%
#     mutate(minimum = min(med_jacc),
#            med = median(med_jacc),
#            maximum = max(med_jacc))

# Robustness Fig 15


## V. House_n ----
house_tibble <- comp_summ_plot(pattern = "house_n",
                               comparison_out = comparison_out,
                               box_name = "SuppFig_20_housekeeping_n.pdf",
                               heat_name = "SuppFig21_house_n_JI_heat.pdf")
house_tibble %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))



## Compile by celltype heatmap plots ----
compile_celltype_heatmaps(pattern="comp_n",
                          relevant_files = "cp_frequencies.RDS",
                          main_title = "Relative\nFrequency",
                          fignum = 22)

compile_celltype_heatmaps(pattern="house_n",
                          relevant_files = "cp_frequencies.RDS",
                          main_title = "Relative\nFrequency",
                          fignum = 23)

compile_celltype_heatmaps(pattern="comp_n",
                          relevant_files = "cp_strength.RDS",
                          main_title = "Relative\nStrength",
                          fignum = 24)

compile_celltype_heatmaps(pattern="house_n",
                          relevant_files = "cp_strength.RDS",
                          main_title = "Relative\nStrength",
                          fignum = 25)

