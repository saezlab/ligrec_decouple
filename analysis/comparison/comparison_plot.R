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

# Comparison out
comparison_out <- "data/output/comparison_out/"


## I. Comp_n ----
comp_tibble <- comp_summ_plot(pattern = "comp_n",
                              comparison_out = comparison_out,
                              box_name = "Figure4.pdf",
                              heat_name = "SuppFig11_mixed_n_JI_heat.pdf")
comp_tibble %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))
# Note CBMC has ~identical Jaccard for different resources and methods



## II. Comp_frac ----
comp_tibble_frac <- comp_summ_plot(pattern = "comp_frac",
                                   comparison_out = comparison_out,
                                   box_name = "SuppFig_10_Composite_frac.pdf",
                                   heat_name = "SuppFig10_mixed_frac_JI_heat.pdf")
comp_tibble_frac %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))


## III. Specs_n ----
n_tibble <- comp_summ_plot(pattern = "specs_n",
                           comparison_out = comparison_out,
                           box_name = "SuppFig_12_Specificity_n.pdf",
                           heat_name = "SuppFig12_specs_n_JI_heat.pdf")
n_tibble %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))


## IV. Specs_FRAC ----
frac_tibble <- comp_summ_plot(pattern = "specs_frac",
                              comparison_out = comparison_out,
                              box_name = "SuppFig_12_Specificity_frac.pdf",
                              heat_name = "SuppFig12_specs_frac_JI_heat.pdf")
frac_tibble %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))

# Robustness Supp. 13


## V. House_n ----
house_tibble <- comp_summ_plot(pattern = "house_n",
                               comparison_out = comparison_out,
                               box_name = "SuppFig_14_housekeeping_n.pdf",
                               heat_name = "SuppFig14_house_n_JI_heat.pdf")
house_tibble %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))


