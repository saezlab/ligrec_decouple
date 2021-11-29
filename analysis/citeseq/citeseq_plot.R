source("analysis/citeseq/citeseq_src.R")
source("analysis/comparison/comparison_utils.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

require(liana)
require(tidyverse)
require(Seurat)
require(ComplexHeatmap)
require(magrittr)


### CURVES ----
settings_vec <- c("specs_n", "comp_n", "house_n")


map(settings_vec, function(score_mode){

    # Generate Plots
    pr_roc_tibble <- readRDS(
        file.path(
            "data/output/citeseq_out/",
            str_glue("citeseq_aurocs_{score_mode}_independent.RDS")
            )
        )

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
        mutate(method_name = recode_methods(method_name)) %>%
        mutate(dataset = recode_datasets(dataset))


    roc_min <- ifelse(min(aucs$roc_auc) > 0.5, 0.5, min(aucs$roc_auc))
    prc_min <- ifelse(min(aucs$prc_auc) > 0.5, 0.5, min(aucs$prc_auc))


    path <- file.path("figures",
                      str_glue("SuppFigX_citeseq_{score_mode}.RDS"))
    cairo_pdf(path,
              height = 15,
              width = 20,
              family = 'DINPro')
    print(
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
            xlim(roc_min, 1) +
            ylim(prc_min, 1) +
            geom_hline(yintercept = 0.5, colour = "pink",
                       linetype = 2, size = 1.2) +
            geom_vline(xintercept = 0.5, colour = "pink",
                       linetype = 2, size = 1.2) +
            theme_bw(base_size = 30) +
            guides(shape=guide_legend(title="Dataset"),
                   color=guide_legend(title="Method"))
    )
    dev.off()
})


# Generate Plots
pr_roc_tibble <- readRDS("data/output/citeseq_out/citeseq_aurocs_independent.RDS")

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



# # AUROC
# auroc_tib <- get_auroc_heat(pr_roc_tibble,
#                             "roc",
#                             mat_only = TRUE)
# pairwise_contrasts <- ggpubr::compare_means(estimate ~ method,
#                                             data = auroc_tib,
#                                             method = "t.test") %>%
#     filter(p.adj <= 0.05)
# pairwise_contrasts
# pairwise_contrasts %>%
#     dplyr::select(group1, group2, p, p.adj, p.signif) %>%
#     as.data.frame() %>%
#     write.csv("~/Downloads/auroc_specificity.csv", row.names = FALSE)
# get_auroc_heat(pr_roc_tibble, "roc", auc_min = 0, auc_max = 1)
#
# # AUPRC
# auprc_tib <- get_auroc_heat(pr_roc_tibble, "prc", mat_only = TRUE)
# pairwise_contrasts <- ggpubr::compare_means(estimate ~ method,
#                                             data = auprc_tib,
#                                             method = "t.test") %>%
#     filter(p.adj <= 0.05)
# pairwise_contrasts
# pairwise_contrasts %>%
#     dplyr::select(group1, group2, p, p.adj, p.signif) %>%
#     as.data.frame() %>%
#     write.csv("~/Downloads/auprc_specificity.csv", row.names = FALSE)
# get_auroc_heat(pr_roc_tibble, "prc", auc_min = 0, auc_max = 1)




## Correlations plot ----
# Load and plot
corr_table <- readRDS("data/output/citeseq_out/citeseq_correlations.RDS")
corr_table %<>%
    # remove Mean/Median ranks
    filter(!(method %in% c("median_rank", "mean_rank"))) %>%
    # rename methods
    mutate(method = gsub("\\..*","", method)) %>%
    # recode methods
    mutate(method = recode_methods(method)) %>%
    # recode datasets
    mutate(dataset = recode_datasets(dataset)) %>%
    # this thing has too few proteins and is an outlier
    filter(dataset!="3kCBMCs")


pairwise_contrasts <- ggpubr::compare_means(estimate ~ method,
                                            data = corr_table,
                                            method = "t.test")
my_comparisons <- pairwise_contrasts %>%
    filter(p.adj <= 0.05) %>%
    select(group1, group2) %>%
    rowwise() %>%
    mutate(my_comparisons = list(c(group1, group2))) %>%
    pluck("my_comparisons")

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05),
                    symbols = c("****", "***", "**", "*"))

# Bar plot
corr_table %>%
    ggplot(aes(x = factor(method), y = estimate)) +
    geom_boxplot(aes(fill = method), alpha = 0.15, show.legend = FALSE) +
    geom_jitter(aes(color = method, shape = dataset, size = n_genes))  +
    scale_shape_manual(values = rep(4:12, len = 7)) +
    xlab("") +
    ylab("Kendal's tau Coefficient") +
    theme_minimal(base_size = 26) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 28)) +
    guides(size=guide_legend(title="Receptor genes"),
           fill= "none",
           color = "none") +
    labs(fill=guide_legend(title="Method"),
         shape=guide_legend(title="Dataset")) +
    geom_signif(comparisons = my_comparisons,
                map_signif_level=TRUE)

