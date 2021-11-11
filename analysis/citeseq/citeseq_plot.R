source("analysis/citeseq/citeseq_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

require(liana)
require(tidyverse)
require(Seurat)
require(ComplexHeatmap)
require(magrittr)

# Generate Plots
pr_roc_tibble <- readRDS("data/output/citeseq_out/citeseq_aurocs_intersect.RDS")

# AUROC
auroc_tib <- get_auroc_heat(pr_roc_tibble,
                            "roc",
                            mat_only = TRUE)
pairwise_contrasts <- ggpubr::compare_means(estimate ~ method,
                                            data = auroc_tib,
                                            method = "t.test") %>%
    filter(p.adj <=0.05)
pairwise_contrasts
pairwise_contrasts %>%
    dplyr::select(group1, group2, p, p.adj, p.signif) %>%
    as.data.frame() %>%
    write.csv("~/Downloads/auroc_specificity.csv", row.names = FALSE)
get_auroc_heat(pr_roc_tibble, "roc", auc_min = 0, auc_max = 1)

# AUPRC
auprc_tib <- get_auroc_heat(pr_roc_tibble, "prc", mat_only = TRUE)
pairwise_contrasts <- ggpubr::compare_means(estimate ~ method,
                                            data = auprc_tib,
                                            method = "t.test") %>%
    filter(p.adj <= 0.05)
pairwise_contrasts
pairwise_contrasts %>%
    dplyr::select(group1, group2, p, p.adj, p.signif) %>%
    as.data.frame() %>%
    write.csv("~/Downloads/auprc_specificity.csv", row.names = FALSE)
get_auroc_heat(pr_roc_tibble, "prc", auc_min = 0, auc_max = 1)


