require(tidyverse)
require(Seurat)
require(liana)

liana_res <- readRDS("data/output/crc_all_test.RDS")
liana_agg <- liana_res %>% liana_aggregate()
liana_agg
