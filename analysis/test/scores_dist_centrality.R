require(liana)
require(tidyverse)
require(magrittr)
require(RColorBrewer)
require(pheatmap)
require(proxy)
require(UpSetR)

source("analysis/comparison/comparison_utils.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

# LOAD New Results from different Method-Resource Combinations
liana_all <- readRDS("data/output/liana_all_resources.RDS")


liana_all_spec <- get_spec_list("data/output/liana_all_resources.RDS",
    #"data/output/crc_res/liana_crc_all.rds", # crc results
                                .score_spec = liana:::.score_specs) # liana:::.score_*

# Obtain Fractions (here we obtain all interactions)

# Distributions of full methods - no filtering (only DE for Connectome)
# As done in the comparisons - i.e. the one to present
plot_score_distributions(liana_all_spec,
                         hit_prop = 1,
                         resource = "OmniPath",
                         pval_thresh = 1,
                         sca_thresh = 0,
                         de_thresh = 0.05)

# PCA by ranks


