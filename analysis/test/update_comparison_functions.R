require(liana)
require(tidyverse)
require(magrittr)
require(RColorBrewer)
require(pheatmap)

source("src/comparison_utils.R")


# Generate New Output
liana_path <- system.file(package = "liana")
testdata <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

liana_all <- liana_wrap(testdata,
                        resource = liana::show_resources()[-c(1:13)], # [-c(1:2)] all resources except Default
                        expr_prop = 0)
saveRDS(liana_all, "data/output/liana_all_resources.RDS")


# LOAD New Results from different Method-Resource Combinations
liana_all <- readRDS("data/output/liana_all_resources.RDS")





#### Load OLD Results from different Method-Resource Combinations
# We only define the measures from each method that we wish to have in our analysis
spec_list <- list("CellChat" =
                      methods::new("MethodSpecifics",
                                   method_name="CellChat",
                                   method_results = readRDS("data/output/crc_res/cellchat_results.rds"),
                                   method_scores=list(
                                       "prob"=TRUE
                                   )),
                  "Connectome" =
                      methods::new("MethodSpecifics",
                                   method_name="Connectome",
                                   method_results = readRDS("data/output/crc_res/conn_results.rds"),
                                   method_scores=list(
                                       "weight_sc"=TRUE
                                   )),
                  "iTALK" =
                      methods::new("MethodSpecifics",
                                   method_name="iTALK",
                                   method_results = readRDS("data/output/crc_res/italk_results.rds"),
                                   method_scores=list(
                                       "weight_comb"=TRUE
                                   )),
                  "NATMI" =
                      methods::new("MethodSpecifics",
                                   method_name="NATMI",
                                   method_results = readRDS("data/output/crc_res/natmi_results.rds"),
                                   method_scores=list(
                                       "edge_specificity"=TRUE
                                   )),
                  "SCA" = methods::new("MethodSpecifics",
                                       method_name="SCA",
                                       method_results = readRDS("data/output/crc_res/sca_results.rds"),
                                       method_scores=list(
                                           "LRscore"=TRUE
                                       )),
                  "Squidpy" =
                      methods::new("MethodSpecifics",
                                   method_name="Squidpy",
                                   method_results = readRDS("data/output/crc_res/squidpy_results.rds"),
                                   method_scores=list(
                                       "pvalue"=FALSE
                                   )))


