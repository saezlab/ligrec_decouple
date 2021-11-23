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

## I. Specs_n ----
comp_summ_plot(pattern = "specs_n",
               comparison_out = comparison_out,
               box_name = "Figure4.pdf")

## II. Specs_FRAC ----
comp_summ_plot(pattern = "specs_frac",
               comparison_out = comparison_out,
               box_name = "SuppFig_10_Specificity_frac_complex.pdf")

## III. House_FRAC ----
comp_summ_plot(pattern = "house_frac", # RENAME
               comparison_out = comparison_out,
               box_name = "SuppFig_10_Specificity_housekeeping_frac_complex.pdf")

