require(tidyverse)
require(magrittr)
require(liana)
source("src/eval_utils.R")


# I) Aggragate helper test `liana_aggregate_enh`
# Generate input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
liana_res <- liana_wrap(seurat_object,
                        method=c("squidpy", "cellchat", "call_sca",
                                 "call_natmi",  "logfc", "call_connectome"))
saveRDS(liana_res, "data/default_liana.rds")

# Read and aggragate
liana_res <- readRDS("data/default_liana.rds")
liana_res$logfc <- NULL
liana_house <- liana_res %>% liana_aggregate_enh(.score_mode=liana:::.score_housekeep)
liana_spec <- liana_res %>% liana_aggregate_enh()


# II) Defaults to be used throughout:
def_new <- liana_wrap(seurat_object,
                      method = c('call_natmi', 'call_connectome', 'logfc',
                                 'cellchat', 'call_sca', 'squidpy', "cytotalk"),
                      # this is passed only to squidpy, cellchat, cytotalk, and logfc
                      expr_prop=0.1,
                      cellchat.params=list(nboot=1000,
                                           expr_prop = 0)) # by default as in CellChat
# CHANGES ONLY FOR MURINE RESOURCE!
# II) Defaults to be used throughout:
def_new <- liana_wrap(seurat_object,
                      resource = "custom",
                      external_resource = murine_resource,
                      method = c('call_natmi', 'call_connectome', 'logfc',
                                 'cellchat', 'call_sca', 'squidpy'),
                      expr_prop=0,
                      squidpy.params=list(threshold = 0.1),
                      cellchat.params=list(organism="mouse", nboot=1000))






