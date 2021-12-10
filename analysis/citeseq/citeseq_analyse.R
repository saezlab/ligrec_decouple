source("analysis/citeseq/citeseq_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")
source("analysis/comparison/comparison_utils.R")

require(liana)
require(tidyverse)
require(magrittr)
require(Seurat)
require(ComplexHeatmap)
require(ggsignif)
require(yardstick)
require(SingleCellExperiment)

# load prereqs
citeseq_dir <- "data/input/citeseq/"
op_resource <- select_resource("OmniPath")[[1]]
murine_resource <- readRDS("data/input/murine_omnipath.RDS")
arbitrary_thresh = 1.645 # one-tailed alpha = 0.05

# Unchanged variables:
.eval = c("max", "independent")
correlation <- TRUE

# Different settings to use
.setting <- c("specs_n",
             "comp_n"#,
             #"house_n"
             )  # n makes no difference
combinations <- expand_grid(.eval, .setting)
combinations


pmap(combinations, function(.eval, .setting){
    print(paste(.eval, .setting))
    set_aggregation_settings(.setting)

    ### Receptor Specificity ROC -----
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
                             adt_pipe_type = "specificity",
                             # liana_aggregate_enh params
                             filt_de_pvals = TRUE,
                             de_thresh = de_thresh, # we only filter Connectome DEs
                             filt_outs = TRUE,
                             pval_thresh = pval_thresh,
                             sca_thresh = 0,
                             .score_mode = .score_specs,
                             .eval = .eval
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
                             adt_pipe_type = "specificity",
                             # liana_aggregate_enh params
                             filt_de_pvals = TRUE,
                             de_thresh = de_thresh, # we only filter Connectome DEs
                             filt_outs = TRUE,
                             pval_thresh = pval_thresh,
                             sca_thresh = 0,
                             .score_mode = .score_specs,
                             .eval = .eval

                )
            }
        }) %>%
        setNames(list.files(citeseq_dir)) %>%
        enframe(name = "dataset") %>%
        unnest(value)

    # Save obj
    saveRDS(pr_roc_tibble, str_glue("data/output/citeseq_out/citeseq_aurocs_{.setting}_{.eval}.RDS"))


    ### Supp) Correlations ----
    if(correlation){
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
                                 adt_pipe_type = "correlation",
                                 # liana_aggregate_enh params
                                 filt_de_pvals = TRUE,
                                 de_thresh = de_thresh, # we only filter Connectome DEs
                                 filt_outs = TRUE,
                                 pval_thresh = pval_thresh,
                                 sca_thresh = 0,
                                 .score_mode = .score_specs
                    )
                } else { # human
                    run_adt_pipe(dir = citeseq_dir,
                                 subdir = subdir,
                                 op_resource = op_resource,
                                 cluster_key = "seurat_clusters",
                                 liana_pattern = "liana_res-0.1.RDS",
                                 organism = "human",
                                 adt_pipe_type = "correlation",
                                 # liana_aggregate_enh params
                                 filt_de_pvals = TRUE,
                                 de_thresh = de_thresh, # we only filter Connectome DEs
                                 filt_outs = TRUE,
                                 pval_thresh = pval_thresh,
                                 sca_thresh = 0,
                                 .score_mode = .score_specs
                    )
                }

            }) %>%
            setNames(list.files(citeseq_dir)) %>%
            enframe(name = "dataset") %>%
            unnest(value)
        saveRDS(corr_table, str_glue("data/output/citeseq_out/citeseq_correlations_{.setting}_{.eval}.RDS"))
        }
    })
