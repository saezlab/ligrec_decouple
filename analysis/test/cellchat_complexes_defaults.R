# issue with CellChat complex
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

# Generate New Output
liana_path <- system.file(package = "liana")
testdata <- readRDS(file.path(liana_path , "testdata",
                              "input", "testdata.rds"))

# Fix Default
liana_all <- liana_wrap(testdata,
                        resource = liana::show_resources()[c(2:6)], # [-c(1:2)] all resources except Default
                        method = c('call_natmi', 'call_connectome', 'logfc',
                                   'cellchat', 'call_sca', 'squidpy'),
                        expr_prop=0)
saveRDS(liana_all, "data/output/liana_default_resources.RDS")
