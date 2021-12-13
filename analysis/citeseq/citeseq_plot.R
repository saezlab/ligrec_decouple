source("analysis/citeseq/citeseq_src.R")
source("analysis/comparison/comparison_utils.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

require(liana)
require(tidyverse)
require(Seurat)
require(ComplexHeatmap)
require(magrittr)
require(patchwork)

### CURVES ----
score_mode <- c("comp_n", "specs_n")
.eval <- c("independent", "intersect", "max")

plot_tib <- expand_grid(score_mode, .eval) %>%
    arrange(.eval)
plot_tib

plot_tib %<>% pmap(function(score_mode, .eval){

    # Generate Plots
    pr_roc_tibble <- readRDS(
        file.path(
            "data/output/citeseq_out/",
            str_glue("citeseq_aurocs_{score_mode}_{.eval}.RDS")
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

    p <- ggplot(aucs,
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

    pname <- str_glue("citeseq_{score_mode}_{.eval}")

    # Bind together
    tibble_row(pname, p)

    }) %>%
    bind_rows()


# Assemble and Print
path <- file.path("figures",
                  "SuppFig24_citeseq.pdf")
pp <- patchwork::wrap_plots(
    as.list(plot_tib$p),
    ncol=1,
    nrow(4)) +
    plot_annotation(tag_levels = 'A', tag_suffix = ')') &
    theme(plot.tag = element_text(face = 'bold', size = 42))
cairo_pdf(path,
          height = 65, #45
          width = 18,
          family = 'DINPro')
print(pp)
dev.off()


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

