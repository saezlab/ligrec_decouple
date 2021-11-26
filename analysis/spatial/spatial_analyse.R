# Script to prep data for main spatial Figures
# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)

source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")



# Loop over all combinations
# eval does not matter here, unless it's intersect
eval_vec <- c("independent"#,
              #"max",
              #"intersect"
              )
score_mode_vec <- c("mixed",
                    # "house",
                    "specs"
                    )

comb_tibble <- expand_grid(eval_vec, score_mode_vec)
.eval = "independent"
score_mode = "specs"


# MOUSE BRAIN ATLAS ----
brain_dir <- "data/input/spatial/brain_cortex/"
# murine_resource <- readRDS("data/input/murine_omnipath.RDS")

# Load Liana
liana_format <- readRDS(str_glue("data/output/aggregates/brain_{.eval}_{score_mode}.RDS")) %>%
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
        decon_mtrx <- deconv_res[[2]]
        # correlations of proportions
        decon_cor <- cor(decon_mtrx)

        # format and z-transform deconv proportion correlations
        deconv_corr_long <- decon_cor %>%
            reshape_coloc_estimate(z_scale = FALSE)

        return(deconv_corr_long)
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


############################# Breast Cancer ####################################
brca_dir <- "data/input/spatial/Wu_etal_2021_BRCA"
visium_dict <- list("1142243F" = "TNBC",
                    "1160920F" = "TNBC",
                    "CID4290" = "ER",
                    "CID4465" = "TNBC",
                    "CID4535" = "ER",
                    "CID44971" = "TNBC")

# tibble over which we pmap
tibble_dict <- tibble(slide_name = names(visium_dict),
                      slide_subtype = visium_dict) %>%
    unnest(slide_subtype) %>%
    mutate(cluster_key = "celltype_major")
tibble_dict %<>%
    bind_rows(tibble_dict %>%
                  mutate(cluster_key = "celltype_minor"))


# NOTE WE KEEP ONLY celltype minor (major are also random) :D
tibble_dict %<>% filter(cluster_key == "celltype_minor")

# Get Deconvolution results
deconv_results <- tibble_dict %>%
    mutate(deconv_results = pmap(tibble_dict, function(slide_name,
                                                       slide_subtype,
                                                       cluster_key){
        # deconvolution subdirectory
        deconv_directory <- file.path(brca_dir,
                                      "deconv", str_glue("{slide_subtype}_{cluster_key}"))

        # load deconv results
        deconv_res <- readRDS(file.path(deconv_directory,
                                        str_glue("{slide_name}_{slide_subtype}_{cluster_key}_deconv.RDS")))
        decon_mtrx <- deconv_res[[2]]

        # correlations of proportions
        decon_cor <- cor(decon_mtrx)

        # format and z-transform deconv proportion correlations
        deconv_corr_long <- decon_cor %>%
            reshape_coloc_estimate(z_scale = FALSE) %>%
            mutate(estimate = replace_na(estimate, 0))

        return(deconv_corr_long)
    }))



# load LIANA results and get coloc
lr_coloc <- pmap(.l = list(deconv_results$slide_name,
                           deconv_results$slide_subtype,
                           deconv_results$deconv_results # deconvolution results
                           ),
           .f = function(slide_name,
                         slide_subtype,
                         deconv){

               message(slide_subtype)

               # load aggregated liana
               liana_format <-
                   readRDS(
                       file.path(
                       "data/output/aggregates/",
                       str_glue("{slide_subtype}_{.eval}_{score_mode}.RDS")
                       )) %>%
                   liana_agg_to_long()
               gc()

               # Format further
               liana_loc <- liana_format %>%
                   left_join(deconv, by = c("source"="celltype1",
                                            "target"="celltype2")) %>%
                   # FILTER AUTOCRINE
                   filter(source!=target) %>%
                   ungroup() %>%
                   mutate(dataset = slide_name) %>%
                   mutate(subtype = slide_subtype)
           }) %>%
    bind_rows()
saveRDS(lr_coloc, str_glue("data/output/spatial_out/Wu_etal_2021_BRCA/coloc_{.eval}_{score_mode}.RDS"))
