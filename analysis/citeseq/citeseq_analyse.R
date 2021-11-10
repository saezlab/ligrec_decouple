source("analysis/citeseq/citeseq_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

require(liana)
require(tidyverse)
require(magrittr)
require(Seurat)
require(ComplexHeatmap)
require(ggsignif)
require(yardstick)
require(SingleCellExperiment)
arbitrary_thresh = 1.645 # one-tailed alpha = 0.05

# load prereqs
citeseq_dir <- "data/input/citeseq/"
op_resource <- select_resource("OmniPath")[[1]]
murine_resource <- readRDS("data/input/murine_omnipath.RDS")

### I) Correlations ----
corr_table <- list.files(citeseq_dir) %>%
    map(function(subdir){
        # If mouse, load convert use murine-specific conversion
        if(stringr::str_detect(subdir, pattern = "spleen")){
            run_adt_pipe(dir = citeseq_dir,
                         subdir = subdir,
                         op_resource = murine_resource,
                         cluster_key = "seurat_clusters",
                         liana_pattern = "liana_res-0.1.RDS",
                         organism = "mouse",
                         adt_pipe_type = "correlation"
            )
        } else { # human
            run_adt_pipe(dir = citeseq_dir,
                         subdir = subdir,
                         op_resource = op_resource,
                         cluster_key = "seurat_clusters",
                         liana_pattern = "liana_res-0.1.RDS",
                         organism = "human",
                         adt_pipe_type = "correlation"
            )
        }

    }) %>% setNames(list.files(citeseq_dir)) %>%
    enframe() %>%
    unnest(value) %>%
    rename(dataset = name)
saveRDS(corr_table, "data/output/citeseq_out/citeseq_correlations.RDS")

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



### II) Receptor Specificity ROC -----
pr_roc_tibble <- list.files(citeseq_dir) %>%
    map(function(subdir){
        # If mouse, load convert use murine-specific conversion
        if(stringr::str_detect(subdir, pattern = "spleen")){
            run_adt_pipe(subdir = subdir,
                         dir = citeseq_dir,
                         op_resource = murine_resource,
                         organism = "mouse",
                         cluster_key = "seurat_clusters",
                         sobj_pattern = "_seurat.RDS",
                         liana_pattern = "liana_res-0.1.RDS",
                         arbitrary_thresh = arbitrary_thresh,
                         adt_pipe_type = "specificity"
            )
        } else { # human
            run_adt_pipe(subdir = subdir,
                         dir = citeseq_dir,
                         op_resource = op_resource,
                         organism = "human",
                         cluster_key = "seurat_clusters",
                         sobj_pattern = "_seurat.RDS",
                         liana_pattern = "liana_res-0.1.RDS",
                         arbitrary_thresh = arbitrary_thresh,
                         adt_pipe_type = "specificity"
            )
        }

    }) %>% setNames(list.files(citeseq_dir)) %>%
    enframe(name = "dataset") %>%
    unnest(value)

# Save obj
saveRDS(pr_roc_tibble, "data/output/citeseq_out/citeseq_aurocs.RDS")

