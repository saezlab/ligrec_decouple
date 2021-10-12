# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)

source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")
brain_dir <- "data/input/spatial/brain_cortex/"
arb_thresh = 1.645 # one-tailed alpha = 0.05

# Load Liana
liana_res <- readRDS(str_glue("{brain_dir}/brain_liana_results.RDS"))
liana_res %<>% liana_aggregate()

# Format LIANA res to long tibble /w method_name and predictor (ranks)
liana_format <- liana_res %>%
    mutate(across(c(source, target), ~str_replace_all(.x, " ", "."))) %>%
    dplyr::select(source, target, ends_with("rank"), -c(mean_rank, median_rank)) %>%
    mutate(aggregate_rank = min_rank(aggregate_rank)) %>%
    pivot_longer(-c(source,target),
                 names_to = "method_name",
                 values_to = "predictor") #rank


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
            reshape_coloc_estimate(z_scale = TRUE)

        return(deconv_corr_long)
    }) %>%
    setNames(slides)




# I) Enrichment of interactions between co-localized cells ----
n_ranks = c(100, 250, 500, 1000, 5000, 10000, 50000)



# map over deconvolution correlations for each slide
fets <- deconv_results %>%
    map2(names(.), function(deconv_cor_formatted, slide_name){
        # Assign colocalisation to liana results in long according to a threshold
        liana_loc <- liana_format %>%
            left_join(deconv_cor_formatted, by = c("source"="celltype1",
                                                   "target"="celltype2"))  %>%
            # FILTER AUTOCRINE
            filter(source!=target) %>%
            dplyr::mutate(localisation = case_when(estimate >= arb_thresh ~ "colocalized",
                                                   estimate < arb_thresh ~ "not_colocalized"
            )) %>% ungroup()

        fets <- n_ranks %>%
            map(function(n_rank){
                run_coloc_fet(liana_loc = liana_loc,
                              n_rank = n_rank) %>%
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

# plot Enrichment of colocalized in top vs total
ggplot(boxplot_data,
            aes(x = n_rank, y = enrichment,
                color = method_name)) +
    geom_boxplot(alpha = 0.15,
                 outlier.size = 1.5,
                 width = 0.2,
                 show.legend = FALSE)  +
    geom_jitter(aes(shape = dataset), width = 0) +
    facet_grid(~method_name, scales='free_x', space='free', switch="x") +
    theme_bw(base_size = 20) +
    geom_hline(yintercept = 0, colour = "red",
               linetype = 2, size = 0.9) +
    theme(strip.text.x = element_text(angle = 90),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
          ) +
    labs(shape=guide_legend(title="Dataset")) +
    ylab("Enrichment") +
    xlab("Number of Ranks Considered")


# II) AUROC/AUPRC Curves -----
# Positive Correlation AUC ----
rocs_positive <- deconv_results %>%
    map2(names(.), function(deconv_cor_formatted, slide_name){
        # Assign colocalisation to liana results in long according to a threshold
        liana_loc <- liana_format %>%
            left_join(deconv_cor_formatted, by = c("source"="celltype1",
                                                   "target"="celltype2"))  %>%
            # FILTER AUTOCRINE
            filter(source!=target) %>%
            # Assign Positive and Negative Classes according to a threshold
            dplyr::mutate(response = case_when(estimate >= arb_thresh ~ 1,
                                               estimate < arb_thresh ~ 0)) %>%
            unite(source, target, col = "celltype_pair") %>%
            mutate(response = as.factor(response)) %>%
            mutate(predictor = predictor * -1)

        auc_df <- liana_loc %>%
            group_by(method_name) %>%
            group_nest(.key = "method_res") %>%
            mutate(prc = map(method_res,
                             function(df)
                                 calc_curve(df,
                                            curve = "PR",
                                            downsampling = TRUE,
                                            times = 100,
                                            source_name = "interaction"))) %>%
            mutate(roc = map(method_res,
                             function(df) calc_curve(df,
                                                     curve="ROC",
                                                     source_name = "interaction"))) %>%
            mutate(dataset = slide_name)

        # prevent RAM from exploding...
        gc()

        return(auc_df)
        }) %>%
    bind_rows()
saveRDS(rocs_positive, "data/output/spatial_out/brain_cortex/rocs_positive.RDS")

# Load results
# ROC
rocs_positive <- readRDS("data/output/spatial_out/brain_cortex/rocs_positive.RDS")
pos_roc <- get_auroc_heat(rocs_positive, "roc") # all random

# PRC
pos_prc <- get_auroc_heat(rocs_positive, "prc")
gc()


# Negative Correlation AUC ----
# Assign Positive and Negative Classes according to a threshold
rocs_negative <- deconv_results %>%
    map2(names(.), function(deconv_cor_formatted, slide_name){
        # Assign colocalisation to liana results in long according to a threshold
        liana_loc <- liana_format %>%
            left_join(deconv_cor_formatted, by = c("source"="celltype1",
                                                   "target"="celltype2"))  %>%
            # FILTER AUTOCRINE
            filter(source!=target) %>%
            # Assign Positive and Negative Classes according to a threshold
            dplyr::mutate(response = case_when(estimate <= -arb_thresh ~ 1,
                                               estimate > arb_thresh ~ 0)) %>%
            unite(source, target, col = "celltype_pair") %>%
            mutate(response = as.factor(response)) %>%
            mutate(predictor = predictor)

        auc_df <- liana_loc %>%
            group_by(method_name) %>%
            group_nest(.key = "method_res") %>%
            mutate(prc = map(method_res,
                             function(df)
                                 calc_curve(df,
                                            curve = "PR",
                                            downsampling = TRUE,
                                            times = 100,
                                            source_name = "interaction"))) %>%
            mutate(roc = map(method_res,
                             function(df) calc_curve(df,
                                                     curve="ROC",
                                                     source_name = "interaction"))) %>%
            mutate(dataset = slide_name)

        # prevent RAM from exploding...
        gc()

        return(auc_df)
    }) %>%
    bind_rows()
saveRDS(rocs_negative, "data/output/spatial_out/brain_cortex/rocs_negative.RDS")

# ROC
rocs_negative <- readRDS("data/output/spatial_out/brain_cortex/rocs_negative.RDS")
neg_roc <- get_auroc_heat(rocs_negative, "roc") # all random
neg_roc

# PRC
neg_prc <- get_auroc_heat(rocs_negative, "prc")
neg_prc
