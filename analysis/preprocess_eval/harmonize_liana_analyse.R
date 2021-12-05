# Script to Re-do the Cytosig and Visium-Spatial evaluations
# with harmonized filtering


# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)
require(SingleCellExperiment)
require(decoupleR)

source("analysis/cytosig/cytosig_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

### Cytosig ----
.eval = "independent"
# Customize tibble to harmonized results
path_tibble <-
    tibble(dataset = c("ER",
                       "HER2",
                       "TNBC"),
           # BRCA seurat objects used throughout the manuscript
           seurat_path = c("data/input/spatial/Wu_etal_2021_BRCA/deconv/ER_celltype_minor/ER_celltype_minor_seurat.RDS",
                           "data/input/spatial/Wu_etal_2021_BRCA/deconv/HER2_celltype_minor/HER2_celltype_minor_seurat.RDS",
                           "data/input/spatial/Wu_etal_2021_BRCA/deconv/TNBC_celltype_minor/TNBC_celltype_minor_seurat.RDS"
           ),
           # liana results generated from extract_evals.sh
           liana_path = c(file.path("data/output/eval_harmonize", str_glue("BRCA_ER_liana_agg.RDS")),
                          file.path("data/output/eval_harmonize", str_glue("BRCA_HER2_liana_agg.RDS")),
                          file.path("data/output/eval_harmonize", str_glue("BRCA_TNBC_liana_agg.RDS"))
           ))

# Run Cytosig pipe
cytosig_eval_wrap(.eval = .eval,
                  score_mode = "comp",
                  generate = FALSE,
                  path_tibble = path_tibble,
                  outpath = "data/output/eval_harmonize")


# Plot results
plot_cytosig_aucs(inputpath = "data/output/eval_harmonize/cytosig_res_independent_comp.RDS")



### Spatial ----
.eval = "independent"
score_mode = "mixed"

# MOUSE BRAIN ATLAS ----
brain_dir <- "data/input/spatial/brain_cortex/"
# murine_resource <- readRDS("data/input/murine_omnipath.RDS")

# Load Liana
liana_format <- readRDS(str_glue("data/output/aggregates/brain_{.eval}_{score_mode}_liana_res.RDS")) %>%
    liana_agg_to_long()

# load deconvolution results and do correlation
slides <- c("anterior1",
            "anterior2",
            "posterior1",
            "posterior2")
deconv_results <- slides %>%
    map(function(slide){
        # load results
        deconv_res <- readRDS(str_glue("{brain_dir}/{slide}_doconvolution2.RDS"))
        deconv_res[[2]]
        # # correlations of proportions
        # decon_cor <- cor(decon_mtrx)
        #
        # # format and z-transform deconv proportion correlations
        # deconv_corr_long <- decon_cor %>%
        #     reshape_coloc_estimate(z_scale = FALSE)
        #
        # return(deconv_corr_long)
    }) %>%
    setNames(slides)

# Bind lr and coloc
lr_coloc <- deconv_results %>%
    map2(names(.), function(deconv_cor_formatted, slide_name){
        # Assign colocalisation to liana results in long according to a threshold
        liana_loc <- liana_format %>%
            left_join(deconv_cor_formatted, by = c("source"="celltype1",
                                                   "target"="celltype2")) %>%
            # FILTER AUTOCRINE
            filter(source!=target) %>%
            ungroup() %>%
            mutate(dataset = slide_name)
    }) %>%
    bind_rows()
# save liana LR score-colocalizations
saveRDS(lr_coloc, str_glue("data/output/spatial_out/brain_cortex/coloc_{.eval}_{score_mode}.RDS"))
