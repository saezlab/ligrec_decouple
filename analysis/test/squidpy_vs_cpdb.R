source("analysis/citeseq/citeseq_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

require(liana)
require(tidyverse)
require(Seurat)
require(ComplexHeatmap)
require(magrittr)


pbmc <- readRDS("data/input/citeseq/5k_pbmcs/5k_pbmcs_seurat.RDS")
pbmc@meta.data %<>% mutate(seurat_clusters = as.factor(seurat_clusters))
Idents(pbmc) <- pbmc@meta.data$seurat_clusters
liana_res <- liana_wrap(pbmc, method = c("squidpy", "cellphonedb"))

# LIANA
liana_cpdb <- liana_res$cellphonedb %>% na.omit()
squidpy_cpdb <- liana_res$squidpy %>% na.omit()

liana_inter
squidpy_inter

# LIANA
liana_inter <- liana_cpdb %>%
    filter(pvalue <= 0.05) %>%
    select(source, ligand, target, receptor)

squidpy_inter <- squidpy_cpdb  %>%
    filter(pvalue <= 0.05) %>%
    select(source, ligand, target, receptor)
setdiff(liana_inter, squidpy_inter)
