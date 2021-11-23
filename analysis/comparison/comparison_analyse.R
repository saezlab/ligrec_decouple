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

# Comparison out
comparison_out <- "data/output/comparison_out/"


## I. Specs_n ----
# Load Jaccard Index Across Resources using the same Method
dirs <- list.files(comparison_out,
                   pattern="specs_n")

across_resource <- map(dirs, function(d){
    readRDS(file.path(comparison_out, d, "across_resource_ji.RDS"))
    }) %>%
    setNames(dirs) %>%
    enframe(name = "dataset_setting",
            value = "ji_stats") %>%
    unnest(ji_stats) %>%
    mutate(dataset_setting =
               factor(dataset_setting,
                      levels = c("cbmc_specs_n",
                                 "panc8_specs_n",
                                 # "crc_specs_n",
                                 "er_specs_n",
                                 "tnbc_specs_n",
                                 "her2_specs_n")))

# Load Jaccard Index Across Methods using the same Resource
across_method <- map(dirs, function(d){
    readRDS(file.path(comparison_out, d, "across_methods_ji.RDS"))
}) %>%
    setNames(dirs) %>%
    enframe(name = "dataset_setting",
            value = "ji_stats") %>%
    unnest(ji_stats) %>%
    mutate(dataset_setting =
               factor(dataset_setting,
                      levels = c("cbmc_specs_n",
                                 "panc8_specs_n",
                                 # "crc_specs_n",
                                 "er_specs_n",
                                 "tnbc_specs_n",
                                 "her2_specs_n")))



# Plot Jaccard Boxplots
jacc_all_boxplot(across_resource,
                 entity="method")

jacc_all_boxplot(across_method,
                 entity="resource")


## II. Specs_FRAC ----
# Load Jaccard Index Across Resources using the same Method
dirs <- list.files(comparison_out,
                   pattern="specs_frac")

across_resource <- map(dirs, function(d){
    readRDS(file.path(comparison_out, d, "across_resource_ji.RDS"))
}) %>%
    setNames(dirs) %>%
    enframe(name = "dataset_setting",
            value = "ji_stats") %>%
    unnest(ji_stats) %>%
    mutate(dataset_setting =
               factor(dataset_setting,
                      levels = c("cbmc_specs_frac",
                                 "panc8_specs_frac",
                                 # "crc_specs_frac",
                                 "er_specs_frac",
                                 "tnbc_specs_frac",
                                 "her2_specs_frac")))

# Load Jaccard Index Across Methods using the same Resource
across_method <- map(dirs, function(d){
    readRDS(file.path(comparison_out, d, "across_methods_ji.RDS"))
}) %>%
    setNames(dirs) %>%
    enframe(name = "dataset_setting",
            value = "ji_stats") %>%
    unnest(ji_stats) %>%
    mutate(dataset_setting =
               factor(dataset_setting,
                      levels = c("cbmc_specs_frac",
                                 "panc8_specs_frac",
                                 # "crc_specs_frac",
                                 "er_specs_frac",
                                 "tnbc_specs_frac",
                                 "her2_specs_frac")))



# Plot Jaccard Boxplots
jacc_all_boxplot(across_resource,
                 entity="method")

jacc_all_boxplot(across_method,
                 entity="resource")

