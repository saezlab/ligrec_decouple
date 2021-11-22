require(tidyverse)
require(magrittr)
require(liana)
require(RColorBrewer)
require(pheatmap)
require(proxy)
require(UpSetR)
require(grid)
require(ComplexHeatmap)
require(patchwork)
require(ggplotify)

source("analysis/comparison/comparison_utils.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

# Runs all things
# saves plot data
# patchwork

# top 500,000 is the default for liana_aggregate_enh
# consider only looking at some resources - no need for all 20

xxx <- comparison_pipe(input_filepath = "data/output/temp/liana_all_resources.RDS",
                       output_filepath = "test123",
                       resource = "OmniPath",
                       top_x = 0.05,
                       top_fun = "top_frac",
                       .score_specs = liana:::.score_specs,
                       cap_value_str = 1,
                       cap_value_freq = 1,
                       pval_thresh = 1,
                       sca_thresh = 0,
                       de_thresh = 1)

# Pre-defined resource list :) (max 8)
# e.g. > 2019 or something


input_filepath = "data/output/comparison_out/panc8_liana_res.RDS"
output_filepath = "panc8_comp"
resource = "ICELLNET"
top_x = 0.05
top_fun = "top_frac"
.score_specs = liana:::.score_specs
cap_value_str = 1
cap_value_freq = 1
pval_thresh = 1
sca_thresh = 0
de_thresh = 0.05




top_hits_key <- str_glue({"top_{top_x}"})
outpath <- str_glue("data/output/comparison_out/{output_filepath}")
message(str_glue("Creating and Saving in : {outpath}"))
dir.create(outpath, showWarnings=FALSE)

# Ranked Scores according to a set of criteria (here by specificity if available)
liana_all_spec <- get_spec_list(input_filepath,
                                .score_spec = .score_specs)



# Top X proportion of hits (according to the ranking specs above)
top_lists <- get_top_hits(liana_all_spec,
                          n_ints = top_x,
                          top_fun = top_fun,
                          de_thresh = 0.05)

ct_strength <- get_ct_strength(liana_all_spec)

strength_heat <- get_ct_heatmap(ct_strength,
                                cap_value = cap_value_str)
