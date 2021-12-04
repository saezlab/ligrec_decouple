# We use this script to run the evaluations with harmonized pre-processing
library(tidyverse)
library(liana)
library(Seurat)
library(magrittr)

source("src/eval_utils.R")
brca_dir <- "data/input/spatial/Wu_etal_2021_BRCA"

## LIANA Runs ----
# A) BRCA ----
# Load BRCA object from Spatial Directory
map(c("HER2", "ER", "TNBC"), function(brca_subtype){
    deconv_directory <- file.path(brca_dir,
                                  "deconv", str_glue("{brca_subtype}_celltype_minor"))

    seurat_object <- readRDS(file.path(deconv_directory,
                                       str_glue("{brca_subtype}_celltype_minor_seurat.RDS")))

    # Here, we run LIANA but with multiple resources
    liana_res <- liana_wrap(seurat_object,
                            method = c("natmi", "connectome", "logfc",
                                       "sca", "cellphonedb", "cellchat", "cytotalk"
                                       ),
                            expr_prop=0.1,
                            cellchat.params = list(nboot=1000,
                                                   expr_prop = 0.1),
                            assay = "SCT")

    # save LIANA results
    message("Saving LIANA RESULTS")
    saveRDS(liana_res,
            file.path("data/output/eval_harmonize",
                      str_glue("BRCA_{brca_subtype}_liana_res.RDS")))
    message("Saving Aggregating LIANA")
    liana_agg <- liana_res %>%
        liana_aggregate_enh(
            filt_outs = FALSE,
            pval_thresh = 1,
            sca_thresh = 0,
            filt_de_pvals = FALSE,
            de_thresh = 1,
            .score_mode = .score_comp,
            .eval = "independent"
        )
    saveRDS(liana_agg,
            file.path("data/output/eval_harmonize",
                      str_glue("BRCA_{brca_subtype}_liana_agg.RDS")))

})



# B) Mouse Brain ----
# brain analysis directory
brain_dir <- "data/input/spatial/brain_cortex"
cortex_sc <- readRDS(file.path(brain_dir, "allen_cortex_prep.rds"))
murine_resource <- readRDS("data/input/murine_omnipath.RDS")
liana_res <- liana_wrap(cortex_sc,
                        method = c("natmi", "connectome", "logfc",
                                   "sca", "cellphonedb", "cellchat", "cytotalk"),
                        resource = "custom",
                        external_resource = murine_resource,
                        expr_prop=0.1,
                        cellchat.params = list(nboot=1000,
                                               organism="mouse",
                                               expr_prop = 0.1),
                        assay = "SCT")
saveRDS(liana_res, "data/output/eval_harmonize/brain_liana_results.RDS")

liana_agg <- liana_res %>%
    liana_aggregate_enh(
        filt_outs = FALSE,
        pval_thresh = 1,
        sca_thresh = 0,
        filt_de_pvals = FALSE,
        de_thresh = 1,
        .score_mode = .score_comp,
        .eval = "independent"
    )
saveRDS(liana_agg, "data/output/eval_harmonize/brain_liana_agg.RDS")
