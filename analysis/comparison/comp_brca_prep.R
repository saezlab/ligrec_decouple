#### Run LIANA with all original methods (with defaults args) + multiple resources
require(tidyverse)
require(magrittr)
require(Seurat)
require(liana)

# Get Args from std in
args <- commandArgs(trailingOnly=TRUE)

# Path to Project
path_to_project <- args[[1]] # ~/Repos/ligrec_decouple/ (on local)
# Get brca_subtpye
brca_subtype <- args[[2]] # "TNBC", "ER", "HER2"
message(str_glue("Now Running: LIANA with {brca_subtype}"))


# Spatial Deconv Directory (i.e atlas directory)
brca_dir <- "data/input/spatial/Wu_etal_2021_BRCA"
deconv_directory <- file.path(path_to_project, brca_dir,
                              "deconv", str_glue("{brca_subtype}_celltype_minor"))

# Load BRCA object from Spatial Directory
seurat_object <- readRDS(file.path(deconv_directory,
                                   str_glue("{brca_subtype}_celltype_minor_seurat.RDS")))


# Here, we run LIANA but with multiple resources
liana_res <- liana_wrap(seurat_object,
                        method = c('call_natmi',
                                   'call_connectome', 'logfc', 'cellchat',
                                   'call_sca', 'cellphonedb', "cytotalk"
                                   ),
                        resource = c("ICELLNET",
                                    "OmniPath",
                                    "CellChatDB"#,
                                    # "CellTalkDB"
                                     ),
                        # this is passed only to squidpy, cellchat, cytotalk, and logfc
                        expr_prop=0.1,
                        cellchat.params = list(nboot=1000,
                                               expr_prop = 0),
                        call_natmi.params = list(
                            expr_file = str_glue("{brca_subtype}_em.csv"),
                            meta_file = str_glue("{brca_subtype}_metadata.csv"),
                            output_dir = str_glue("brca_{brca_subtype}_results"),
                            reso_name = str_glue("{brca_subtype}_placeholder")
                            ),
                        assay = "SCT") # by default as in CellChat

# save LIANA results
message("Saving LIANA RESULTS")
saveRDS(liana_res,
        file.path(path_to_project, "data/output/comparison_out/",
                  str_glue("BRCA_{brca_subtype}_liana_res.RDS")))



# Aggregate and Save OmniPath's out - to be used in Spatial and CytoSig evals
source("src/eval_utils.R")
lr_omni <- liana_res %>%
    transpose() %>%
    pluck("OmniPath") %>%
    liana_aggregate_enh(filt_de_pvals = TRUE,
                        de_thresh = 0.05, # we only filter connectome DEs
                        filt_outs = FALSE,
                        pval_thresh = 1,
                        sca_thresh = 0,
                        .score_mode = liana:::.score_specs,
                        cap = 500000)
saveRDS(lr_omni,
        file.path(path_to_project, "data/output/comparison_out/",
                  str_glue("BRCA_{brca_subtype}_liana_omni.RDS")))
