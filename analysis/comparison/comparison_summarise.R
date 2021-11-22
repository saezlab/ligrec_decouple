## Script to generate the comparison data presented in the Manuscript
# Arg 1: Path to LIANA results (resource-method) as Input
# Arg 2: Name of output folder

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

# Get Args from std in
args <- commandArgs(trailingOnly=TRUE)

# Get Job name
output_filepath <- args[[2]]  # e.g. "panc8_specs_n" (+ job name)
input_filepath <- args[[3]] # e.g. "data/output/comparison_out/panc8_liana_res.RDS"
setting <- args[[4]]

# load appropriate settings
if(setting=="specs_frac"){
    .score_specs = liana:::.score_specs
    top_fun <- "top_frac"
    top_x <- 0.05
} else if(setting=="house_frac"){
    .score_specs = liana:::.score_housekeep
    top_fun <- "top_frac"
    top_x <- 0.05
} else if(setting=="specs_n"){
    .score_specs = liana:::.score_specs
    top_fun <- "top_n"
    top_x <- 1000
} else{
    stop("Setting is wrong!")
}
message(setting)

# Summarize comparisons
comparison_summary(input_filepath = input_filepath,
                   output_filepath = output_filepath,
                   resource = "OmniPath",
                   top_x = top_x,
                   top_fun = top_fun,
                   .score_specs = .score_specs,
                   cap_value_freq = 1,
                   pval_thresh = 1,
                   sca_thresh = 0,
                   de_thresh = 0.05
                   )
