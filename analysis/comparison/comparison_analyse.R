## Script to generate the comparison data presented in the Manuscript
# Arg 1: Path to LIANA results (resource-method) as Input
# Arg 2: Name of output folder

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

# Get Args from std in
args <- commandArgs(trailingOnly=TRUE)

# Get Job name
input_filepath <- args[[2]] # e.g. "data/output/comparison_out/panc8_liana_res.RDS"
output_filepath <- args[[3]] # e.g. panc8_out

# Summarize comparisons
comparison_summary(input_filepath = input_filepath,
                   output_filepath = output_filepath,
                   resource = "OmniPath",
                   top_x = 0.05,
                   top_fun = "top_frac",
                   .score_specs = liana:::.score_specs,
                   cap_value_str = 1,
                   cap_value_freq = 1,
                   pval_thresh = 1,
                   sca_thresh = 0,
                   de_thresh = 0.05
                   )
