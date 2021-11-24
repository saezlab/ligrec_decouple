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
n_tibble <- comp_summ_plot(pattern = "specs_n",
                           comparison_out = comparison_out,
                           box_name = "Figure4.pdf",
                           heat_name = "SuppFig11_JI_heat.pdf")
n_tibble %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))


## PURGATORY ----
pattern="specs_n"
dirs <- list.files(comparison_out,
                   pattern=pattern)

levels <- map_chr(c("cbmc", "panc8", "crc",
                    "er", "tnbc", "her2"),
                  function(ds){
                      paste(ds, pattern, sep = "_")
                  })

# Merge all Binary DFs and Calculate Jaccard Index for each
# * binary DFs are top X hits, full joined across method-resources
all_binary_dfs <- map(dirs, function(d){
    readRDS(file.path(comparison_out, d, "binary_df.RDS"))
}) %>%
    setNames(dirs) %>%
    enframe(name = "dataset_setting",
            value = "binary") %>%
    mutate(jaccard_mat = binary %>% map(function(b){
        get_heatmap_data(b,
                         sim_dist = "simil",
                         method = "Jaccard",
                         diag = TRUE,
                         upper = TRUE)
    }))

# Calculate jaccard index median across datasets
jaccard_mat_mean <- apply(simplify2array(all_binary_dfs$jaccard_mat),
                          1:2,
                          median)

# Plot Merged JI Heat
merged_ji_heat <- get_simdist_heatmap(binary_df = all_binary_dfs$binary[[1]],
                                      simdif_df = jaccard_mat_mean,
                                      sim_dist = "simil",
                                      method = "Jaccard",
                                      diag = TRUE,
                                      upper = TRUE,
                                      cluster_rows = TRUE,
                                      cluster_columns = TRUE)


## II. Specs_FRAC ----
frac_tibble <- comp_summ_plot(pattern = "specs_frac",
                              comparison_out = comparison_out,
                              box_name = "SuppFig_10_Specificity_frac.pdf")
frac_tibble %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))

## III. House_n ----
house_tibble <- comp_summ_plot(pattern = "house_n",
                               comparison_out = comparison_out,
                               box_name = "SuppFig_12_housekeeping_n.pdf")
house_tibble %>%
    group_by(entity) %>%
    mutate(minimum = min(med_jacc),
           med = median(med_jacc),
           maximum = max(med_jacc))
