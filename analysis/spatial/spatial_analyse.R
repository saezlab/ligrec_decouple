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

### Bind All Conditions
visium_dict <- list("1142243F" = "TNBC",
                    "1160920F" = "TNBC",
                    "CID4290" = "ER",
                    "CID4465" = "TNBC",
                    "CID4535" = "ER",
                    "CID44971" = "TNBC",
                    "anterior1" = "brain",
                    "anterior2" = "brain",
                    "posterior1" = "brain",
                    "posterior2" = "brain")


# 1) Get Deconvolution results ----
tibble_dict <- tibble(slide_name = names(visium_dict),
                      slide_condition = visium_dict) %>%
    unnest(slide_condition)


deconv_results <- tibble_dict %>%
    mutate(deconv_results = pmap(tibble_dict, function(slide_name,
                                                       slide_condition){
        if(slide_condition != "brain"){
            # deconvolution subdirectory
            brca_dir <- "data/input/spatial/Wu_etal_2021_BRCA"
            deconv_directory <- file.path(brca_dir,
                                          "deconv", str_glue("{slide_condition}_celltype_minor"))
            slide_path = file.path(deconv_directory,
                                   str_glue("{slide_name}_{slide_condition}_celltype_minor_deconv.RDS"))
        } else if(slide_condition == "brain"){
            brain_dir <- "data/input/spatial/brain_cortex/"
            slide_path <- file.path(
                brain_dir,
                str_glue("{slide_name}_doconvolution2.RDS")
            )

        }
        # load deconv results
        deconv_res <- readRDS(slide_path)
        deconv_res[[2]]
    }))

map(c("ER", "TNBC", "brain"), function(condition){
    # Bind All slides
    deconv_bound <- deconv_results %>%
        filter(slide_condition==condition) %>%
        pluck("deconv_results") %>%
        map(function(deconv) deconv %>% as.data.frame) %>%
        bind_rows() %>%
        select(-res_ss)
    # save bound
    saveRDS(deconv_results,
            file.path("data/output/spatial_out/deconv_summ",
                      str_glue("{condition}_coloc_bound.RDS")))

    # get correlation
    spatial_correlations <- cor(deconv_bound) %>%
        reshape_coloc_estimate(z_scale = TRUE) %>%
        mutate(estimate = replace_na(estimate))
    # save corr of bound
    saveRDS(spatial_correlations,
            file.path("data/output/spatial_out/deconv_summ",
                      str_glue("{condition}_coloc_corr.RDS")))
})


# 2) Get Co-localized Enrichment ----

# Loop over all combinations
# eval does not matter here, unless it's intersect
.eval <- c(#"independent",
           "max" #,
           #"intersect"
           )
score_mode <- c("mixed",
                # "house",
                "specs"
                )

comb_tibble <- expand_grid(.eval, score_mode)
# .eval = "independent"
# score_mode = "mixed"


# Get FET across all conditions
comb_tibble %>%
    pmap(function(.eval, score_mode){

        # Path Tibble -> to be modified
        coloc_tibble <- tibble::tribble(~condition, # condition name
                                        ~liana_agg_path, # path to aggregated liana
                                        ~spatial_corr_path, # path to lr_corr output
#
#                                         "ER+ BRCA",
#                                         file.path(
#                                             "data/output/aggregates/",
#                                             str_glue("ER_{.eval}_{score_mode}_liana_res.RDS")
#                                             ),
#                                         "data/output/spatial_out/deconv_summ/ER_coloc_corr.RDS",

                                        "TNBC",
                                        file.path(
                                            "data/output/aggregates/",
                                            str_glue("TNBC_{.eval}_{score_mode}_liana_res.RDS")
                                        ),
                                        "data/output/spatial_out/deconv_summ/TNBC_coloc_corr.RDS",

                                        "Brain Cortex",
                                        file.path(
                                            "data", "output", "aggregates",
                                            str_glue("brain_{.eval}_{score_mode}_liana_res.RDS")
                                        ),
                                        "data/output/spatial_out/deconv_summ/brain_coloc_corr.RDS"
        )

        type_lr_coloc <- coloc_tibble %>%
            pmap(function(condition, liana_agg_path, spatial_corr_path){
                get_lr_colocalized(liana_agg_path,
                                   spatial_corr_path,
                                   condition = condition,
                                   corr_thresh = 1.645,
                                   n_ranks = c(100, 250, 500, 1000,
                                               2500, 5000, 10000)
                                   )
            }) %>%
            bind_rows()

        saveRDS(type_lr_coloc, str_glue("data/output/spatial_out/all_fets_{score_mode}.RDS"))
    })

