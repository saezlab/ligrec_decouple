# Load required packages
require(tidyverse)
require(OmnipathR)
require(magrittr)
require(proxy)
require(viridis)
require(RCurl)
require(UpSetR)
require(liana) # !!! Keep liana seperate branch when final revision is done.
require(shadowtext)
require(logger)
require(rlang)

# Source required functions
source("src/plot_utils.R")
source("src/resource_descriptive.R")
source("src/result_utils.R")

## Resource descriptive analysis
# Obtain list with CCC Resources
ligrec <- compile_ligrec_descr()

# Run figure pipeline
descriptive_plots(ligrec)
