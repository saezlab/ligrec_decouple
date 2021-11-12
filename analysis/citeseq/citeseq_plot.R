source("analysis/citeseq/citeseq_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

require(liana)
require(tidyverse)
require(Seurat)
require(ComplexHeatmap)
require(magrittr)

# Generate Plots
pr_roc_tibble <- readRDS("data/output/citeseq_out/citeseq_aurocs_independent.RDS")

# AUROC
auroc_tib <- get_auroc_heat(pr_roc_tibble,
                            "roc",
                            mat_only = TRUE)
pairwise_contrasts <- ggpubr::compare_means(estimate ~ method,
                                            data = auroc_tib,
                                            method = "t.test") %>%
    filter(p.adj <= 0.05)
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


### Boxplot
aucs <- pr_roc_tibble %>%
    dplyr::select(-adt_rank) %>%
    unnest(prc) %>%
    dplyr::select(dataset, method_name, roc, prc_auc = auc) %>%
    distinct() %>%
    unnest(roc) %>%
    dplyr::select(dataset, method_name, roc_auc = auc, prc_auc) %>%
    distinct() %>%
    group_by(method_name) %>%
    mutate(roc_mean = mean(roc_auc),
           prc_mean = mean(prc_auc)) %>%
    ungroup() %>%
    mutate(method_name = gsub("\\..*","", method_name)) %>%
    mutate(method_name = recode_methods(method_name))


roc_min <- ifelse(min(aucs$roc_auc) > 0.5, 0.5, min(aucs$roc_auc))
prc_min <- ifelse(min(aucs$prc_auc) > 0.5, 0.5, min(aucs$prc_auc))

# min_lim <- floor(min(c(aucs$roc, aucs$prc)) * 100)/100
# max_lim <- ceiling(max(c(aucs$roc, aucs$prc)) * 100)/100

ggplot(aucs,
       aes(x=roc_mean,
           y=prc_mean,
           color=method_name)) +
    geom_point(shape = 9, size = 12, alpha=1) +
    scale_shape_manual(values = rep(4:12, len = 20)) +
    geom_point(aes(x = roc_auc,
                   prc_auc,
                   shape=dataset),
               size = 6,
               alpha = 0.3) +
    theme(text = element_text(size=16)) +
    xlab('AUROC') +
    ylab('AUPRC') +
    xlim(roc_min, 0.9) +
    ylim(prc_min, 0.9) +
    geom_hline(yintercept = 0.5, colour = "pink",
               linetype = 2, size = 1.2) +
    geom_vline(xintercept = 0.5, colour = "pink",
               linetype = 2, size = 1.2) +
    theme_bw(base_size = 30) +
    guides(shape=guide_legend(title="Dataset"),
           color=guide_legend(title="Method"))

