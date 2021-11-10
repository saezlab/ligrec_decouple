require(tidyverse)
require(magrittr)
require(liana)
require(RColorBrewer)
require(pheatmap)
require(proxy)
require(UpSetR)
require(grid)
require(ComplexHeatmap)

source("analysis/comparison/comparison_utils.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

# Runs all things
# saves plot data
# patchwork

comparison_pipe(input_filepath <- "data/output/temp/liana_all_resources.RDS",
                output_filepath <- "test123",
                resource = "OmniPath",
                top_x = 0.05,
                top_fun = "top_frac",
                .score_specs = liana:::.score_specs,
                cap_value_str = 1,
                cap_value_freq = 1,
                pval_thresh = 1,
                sca_thresh = 0,
                de_thresh = 0.05)

#' Comparison pipe to be run on each individual dataset
#'
#' @param input_filepath input file path
#' @param output_filepath outpufolder (in comparison_out)
#' @param resource used for distributions plot (it's generated only for 1 resource)
#' @param top_x proportion or top x
#' @param top_fun proportion or top x
#' @param .score_specs type of function aggregate function/score ranking
#' @param cap_value_str cap for strength per CT heatmap
#' @param cap_value_freq cap for strength per CT heatmap
#' @inheritDotParams top_enh and liana_aggregate_enh
comparison_pipe <- function(input_filepath,
                            output_filepath,
                            resource = "OmniPath",
                            top_x = 0.05,
                            top_fun = "top_frac",
                            .score_specs = liana:::.score_specs,
                            cap_value_str = 1,
                            cap_value_freq = 1,
                            ...){

    top_hits_key <- str_glue({"top_{top_x}"})
    outpath <- str_glue("data/output/comparison_out/{output_filepath}")
    dir.create(outpath)

    # LOAD New Results from different Method-Resource Combinations
    liana_all <- readRDS(input_filepath)

    # Ranked Scores according to a set of criteria (here by specificy whenever available)
    liana_all_spec <- get_spec_list(input_filepath,
                                    .score_spec = .score_specs)

    # Top X proportion of hits (according to the ranking specs above)
    top_lists <- get_top_hits(liana_all_spec,
                              n_ints = top_x,
                              top_fun = top_fun,
                              ...)


    # I) Score Distributions -----
    # obtain Per method list
    liana_scores <- get_score_distributions(liana_all_spec,
                                            hit_prop = 1,
                                            resource = resource,
                                            ...)
    saveRDS(liana_scores, str_glue("{outpath}/liana_scores.RDS"))

    score_dist_plot <- plot_score_distributions(liana_scores)
    print(score_dist_plot) # print to check at run time :)

    # II) Interaction Relative Strength per Cell Type -----
    ct_strength <- get_ct_strength(liana_all_spec, ...)
    strength_heat <- get_ct_heatmap(ct_strength, cap_value = cap_value_str)

    # III) Interaction Relative Strength per Cell Type (TOP) -----
    ct_frequncies <- get_ct_frequncies(top_lists[[top_hits_key]])
    freq_heat <- get_ct_heatmap(ct_frequncies, cap_value = cap_value_freq)


    # IV) JI Heatmap (TOP) ----
    jacc_heat <- get_simdist_heatmap(top_lists[[top_hits_key]],
                                     sim_dist = "simil",
                                     method = "Jaccard",
                                     diag = TRUE,
                                     upper = TRUE,
                                     cluster_rows = TRUE,
                                     cluster_columns = TRUE)


    # V) JI Stats/Boxplots ----


}
