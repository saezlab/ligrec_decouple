#### Run LIANA with all original methods (with defaults args) + multiple resources
require(tidyverse)
require(magrittr)
require(Seurat)
require(liana)

source("src/eval_utils.R")

map(c("ER", "TNBC", "HER2"), function(subtype,
                                      resource="OmniPath"){
    liana_res <- readRDS(str_glue("data/output/comparison_out/BRCA_{subtype}_liana_res.RDS")) %>%
        transpose() %>%
        pluck(resource) %>%
        liana_aggregate_enh(filt_de_pvals = TRUE,
                            de_thresh = 0.05, # we only filter connectome DEs
                            filt_outs = FALSE,
                            pval_thresh = 1,
                            sca_thresh = 0,
                            .score_mode = liana:::.score_specs)

    saveRDS(liana_res, str_glue("data/output/comparison_out/BRCA_{subtype}_liana_{resource}.RDS"))
})
