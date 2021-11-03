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
    #"data/output/crc_res/liana_crc_all.rds",
                                .score_spec = liana:::.score_specs) # liana:::.score_*


# Obtain top hits lists
top_lists <- get_top_hits(liana_all_spec,
                          n_ints=c(100,
                                   # 250,
                                   500,
                                   1000))

# Obtain Fractions
top_frac_lists <- get_top_hits(liana_all_spec,
                               n_ints=c(1),
                               top_fun = "top_frac")

# Transpose to resource-method
liana_resmet <- top_frac_lists$top_1 %>%
    transpose()


liana_resmet$OmniPath %>%
    map(function(met_res){

    })

liana_resmet$OmniPath$call.natmi

liana:::.score_specs()[["call.natmi"]]


top_frac_lists$top_1$squidpy$OmniPath %>%
    ggplot(aes(x=pvalue)) +
    geom_density()
