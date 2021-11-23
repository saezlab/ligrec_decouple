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
                                 "her2_specs_n"))
           )

# Calculate Median across Methods (Same Method, Different Resources)
across_resource
median_across_resource <- across_resource %>%
    group_by(dataset_setting) %>%
    summarise(med_jacc = median(jacc))
median_across_resource

# Calculate Median Across Resources (Same Resource, Different Methods)
across_method
median_across_method <- across_method %>%
    group_by(dataset_setting) %>%
    summarise(med_jacc = median(jacc))
median_across_method

# Plot Jaccard Boxplots
method_jacc_box <- jacc_all_boxplot(across_resource,
                                    entity="method",
                                    median_across_resource = median_across_resource,
                                    median_across_method = median_across_method
                                    )

resource_jacc_box <- jacc_all_boxplot(across_method,
                                      entity="resource",
                                      median_across_resource = median_across_resource,
                                      median_across_method = median_across_method)


# Jaccard plots assembled with patchwork
cairo_pdf(str_glue("data/output/comparison_out/SuppFig_10_Specificity_n.pdf"),
          width = 24,
          height = 18,
          family = 'DINPro')
print((method_jacc_box /
           resource_jacc_box) +
          plot_layout(guides = 'collect') +
          plot_annotation(tag_levels = 'A',
                          tag_suffix = ')') &
          theme(plot.tag = element_text(face = 'bold',
                                        size = 36)))
dev.off()


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
method_jacc_box <- jacc_all_boxplot(across_resource,
                                    entity="method")

resource_jacc_box <- jacc_all_boxplot(across_method,
                                      entity="resource")


