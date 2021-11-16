# Load required packages
require(tidyverse)
require(OmnipathR)
require(magrittr)
require(proxy)
require(viridis)
require(RCurl)
require(UpSetR)
require(liana)
require(shadowtext)
require(logger)
require(rlang)
library(patchwork)
library(ggplotify)

# Source required functions
source("src/plot_utils.R")
source("analysis/resource_analysis/resource_descriptive.R")

## Resource descriptive analysis
# Obtain list with CCC Resources
ligrec <- compile_ligrec_descr()
saveRDS(ligrec, "data/input/ligrec.RDS")

# Run figure pipeline
ligrec <- readRDS("data/input/ligrec.RDS")
descript <- descriptive_plots(ligrec)
saveRDS(descript, "data/output/descript.RDS")


# Assemble Main Figure





