# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)

source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")
arb_thresh = 1.645 # one-tailed alpha = 0.05

# Load Liana
liana_res <- readRDS("data/input/spatial/brain_cortex/brain_liana_results.RDS")
liana_res %<>% liana_aggregate()
# Format LIANA res to long tibble /w method_name and predictor (ranks)
liana_format <- liana_res %>%
    mutate(across(c(source, target), ~str_replace_all(.x, " ", "."))) %>%
    dplyr::select(source, target, ends_with("rank"), -c(mean_rank, median_rank)) %>%
    mutate(aggregate_rank = min_rank(aggregate_rank)) %>%
    pivot_longer(-c(source,target),
                 names_to = "method_name",
                 values_to = "predictor") #rank

# I) Enrichment of interactions between co-localized cells ----
# load spotlight deconv example
spotlight_ls <- readRDS("data/spotlight_ls.rds")

# Load deconv proportions and estimate correlation
decon_mtrx <- spotlight_ls[[2]]
decon_cor <- cor(decon_mtrx)

# load deconvolution results and do correlation
deconv_results <- slides %>%
    map(function(slide){
        # load results
        deconv_res <- readRDS(str_glue("{brain_dir}/{slide}_doconvolution2.RDS"))
        decon_mtrx <- deconv_res[[2]]
        # correlations of proportions
        decon_cor <- cor(decon_mtrx)

        # format and z-tranform
        deconv_corr_long <- decon_cor %>%
            reshape_coloc_estimate(z_scale = TRUE)

        return(deconv_corr_long)
    }) %>%
    setNames(slides)


# format deconv proportions
decon_corr_long <- decon_cor %>%
    reshape_coloc_estimate(z_scale = TRUE)
n_ranks = c(100, 250, 500, 1000, 5000, 10000, 50000)

# map over deconvolution correlations for each slide
fets <- deconv_results %>%
    map2(names(.), function(deconv_cor_formatted, slide_name){
        fets <- n_ranks %>%
            map(function(n_rank){
                run_coloc_fet(coloc_estimate = deconv_cor_formatted,
                              liana_format = liana_format,
                              n_rank = n_rank,
                              arb_thresh = arb_thresh) %>%
                    mutate(n_rank = n_rank)
            }) %>%
            bind_rows() %>%
            mutate(dataset = slide_name)
    }) %>%
    bind_rows()


# PLOT
# Squidpy and CellChat return always the same number of interactions
# assign rank to those and only show once
cellchat_squidpy_mins <-  liana_format %>%
    filter(predictor <= min(n_ranks)) %>%
    filter(method_name %in% c("squidpy.rank", "cellchat.rank")) %>%
    group_by(method_name) %>%
    summarise(min_rank = n())

boxplot_data <- fets %>%
    # replace squidpy_cellchat ranks
    left_join(cellchat_squidpy_mins) %>%
    mutate(n_rank = ifelse(is.na(min_rank),
                           n_rank,
                           min_rank
                           )) %>%
    select(-min_rank) %>%
    mutate(n_rank = as.factor(n_rank)) %>%
    distinct_at(.vars = c("method_name", "enrichment", "n_rank", "dataset"),
                .keep_all = TRUE) %>%
    # recode methods
    mutate(method_name = gsub("\\..*","", method_name)) %>%
    mutate(method_name = recode_methods(method_name)) %>%
    # recode datasets
    mutate(dataset = recode_datasets(dataset))


ggplot(boxplot_data,
            aes(x = n_rank, y = enrichment,
                color = method_name)) +
    geom_boxplot(aes(fill = method_name),
                 alpha = 0.15,
                 outlier.size = 1.5,
                 width = 0.2)  +
    geom_jitter(aes(shape = dataset), width = 0) +
    facet_grid(~method_name, scales='free_x', space='free', switch="x") +
    theme_bw(base_size = 20) +
    geom_hline(yintercept = 1, colour = "red", linetype = 2, size = 0.9) +
    theme(strip.text.x = element_text(angle = 90),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
          )


