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
                              heat_name = "SuppFig13_mixed_n_JI_heat.pdf")
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
                                   heat_name = "SuppFig14_mixed_frac_JI_heat.pdf")
comp_tibble_frac %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))


## III. Specs_n ----
n_tibble <- comp_summ_plot(pattern = "specs_n",
                           comparison_out = comparison_out,
                           box_name = "SuppFig_15_Specificity_n.pdf",
                           heat_name = "SuppFig16_specs_n_JI_heat.pdf")
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
                               box_name = "SuppFig_17_housekeeping_n.pdf",
                               heat_name = "SuppFig18_house_n_JI_heat.pdf")
house_tibble %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))

## Compile by celltype heatmap plots ----
relevant_dirs <- list.files(comparison_out, pattern = "house_n")
relevant_dirs <- relevant_dirs[order(relevant_dirs[c(1,2,3,4,6,5)])]

cps <- map(relevant_dirs, function(cdir){
    # readRDS(file.path(comparison_out, cdir, "cp_frequencies.RDS")) %>%
        # get_ct_heatmap(main_title="Relative\nFrequency") %>%
    readRDS(file.path(comparison_out, cdir, "cp_strength.RDS")) %>%
        get_ct_heatmap(main_title="Relative\nStrength") %>%
        as.ggplot()
    }) %>%
    setNames(relevant_dirs)



# Pathwork them -> Print Supp Fig
path <- file.path("figures",
                  "SuppFigure19_CompFreq.RDS")
cairo_pdf(path,
          height = 110,
          width = 30,
          family = 'DINPro')
print(patchwork::wrap_plots(cps) +
          plot_layout(guides = 'keep', ncol = 1) +
          plot_annotation(tag_levels = 'A',
                          tag_suffix = ')') &
          theme(plot.tag = element_text(face = 'bold',
                                        size = 40)))
dev.off()

# Output them seperately in temp

