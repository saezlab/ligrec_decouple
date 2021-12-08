### Assemble Robustness Plots
# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(patchwork)

source("src/plot_utils.R")
source("src/eval_utils.R")
source("analysis/robustness/Code/Utilities/User_Outputs_and_Plots.R")

# Subsampling
subsample <- readRDS("analysis/robustness/Outputs/Cluster_Reshuffling/Boxplot_CR_seurat_pbmc_subset_top100_2021-12-06_21-43.rds") %>%
    format_robustness_plot(., descript = "Cells Removed (%)")
subsample

# Cluster reshuffling
reshuffle <- readRDS("analysis/robustness/Outputs/Cluster_Reshuffling/Boxplot_CR_seurat_pbmc_reshuffle_top100_2021-12-06_21-43.rds") %>%
    format_robustness_plot(., descript = "Reshuffled Cell Labels (%)")
reshuffle

# Resource dilution
dilution <- readRDS("analysis/robustness/Outputs/Resource_Dilution/Boxplot_RD_seurat_pbmc_rand_topo_variable_top100_2021-12-06_21-43.rds") %>%
    format_robustness_plot(., descript = "Interactions Replaced (%)")
dilution

# Indiscriminant resource dilution
indilution <- readRDS("analysis/robustness/Outputs/Resource_Dilution/Boxplot_RD_seurat_pbmc_mod_baseline_rand_topo_variable_top100_2021-12-06_21-43.rds") %>%
    format_robustness_plot(., descript = "Indiscriminantly Replaced Interactions (%)")
indilution


### Assemble Robustness Figure
path <- file.path( "figures", "SuppFig15_robustness.pdf")
pp <- patchwork::wrap_plots(
    list(subsample,
         reshuffle,
         dilution,
         indilution),
    ncol=1,
    nrow(4)) +
    plot_annotation(tag_levels = 'A', tag_suffix = ')') &
    theme(plot.tag = element_text(face = 'bold', size = 32))
cairo_pdf(filename = path,
          height = 42,
          width = 18)
print(pp)
dev.off()



