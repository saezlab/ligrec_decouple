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
descript$size_overlap_combined + theme_bw() +
    theme(
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 23),
        axis.title.y = element_text(size = 34),
        # panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.ticks = element_blank(),
        legend.text = element_text(size=21),
        strip.text.x = element_text(size=19),
        legend.key.size = unit(21, 'mm')
    )

