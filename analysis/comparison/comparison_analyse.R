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


## PURGATORY ----

# Load Jaccard Index Across Resources using the same Method
pattern="specs_n"
dirs <- list.files(comparison_out,
                   pattern=pattern)

levels <- map_chr(c("cbmc", "panc8", "crc",
                    "er", "tnbc", "her2"),
                  function(ds){
                      paste(ds, pattern, sep = "_")
                  })

across_resource <- map(dirs, function(d){
    readRDS(file.path(comparison_out, d, ""))
}) %>%
    setNames(dirs) %>%
    enframe(name = "dataset_setting",
            value = "ji_stats") %>%
    unnest(ji_stats) %>%
    mutate(dataset_setting =
               factor(dataset_setting,
                      levels = levels))



## II. Specs_FRAC ----
comp_summ_plot(pattern = "specs_frac",
               comparison_out = comparison_out,
               box_name = "SuppFig_10_Specificity_frac_complex.pdf")

## III. House_FRAC ----
comp_summ_plot(pattern = "house_n", # RENAME
               comparison_out = comparison_out,
               box_name = "SuppFig_10_Specificity_housekeeping_frac_complex.pdf")

