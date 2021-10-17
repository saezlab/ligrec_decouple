# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)

source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")
arb_thresh = 1.645 # one-tailed alpha = 0.05
murine_resource <- readRDS("data/input/murine_omnipath.RDS")
n_ranks = c(50, 100,
            250, 500,
            1000, 5000,
            10000, 50000)


# merFISH genes
# seqfish_obj <- readRDS("data/input/spatial/fishes/seqFISH_seurat.RDS")

# MOUSE BRAIN ATLAS ----
brain_dir <- "data/input/spatial/brain_cortex/"

# Load Liana
liana_res <- readRDS(str_glue("{brain_dir}/brain_liana_results.RDS"))
# Format LIANA res to long tibble /w method_name and predictor (ranks)
liana_format <- liana_res %>%
    liana_aggregate() %>%
    # filter(ligand %in% rownames(seqfish_obj)) %>%
    # filter(receptor %in% rownames(seqfish_obj)) %>%
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
            reshape_coloc_estimate(z_scale = TRUE)

        return(deconv_corr_long)
    }) %>%
    setNames(slides)




# I) Enrichment of interactions between co-localized cells ----
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
    geom_hline(yintercept = 1, colour = "lightblue",
               linetype = 2, size = 0.9) +
    geom_hline(yintercept = -1, colour = "pink",
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




########################### SeqFISH, merFISH ###################################
fish_dir <- "data/input/spatial/fishes"

# Seurat Object paths
sobj_paths <- list.files(fish_dir, pattern = "_seurat") %>%
    map_chr(function(seurat_path){
        file.path(fish_dir, seurat_path)
        })

# NES paths
nes_paths <- list.files(fish_dir, pattern = "_nes") %>%
    map_chr(function(nes_path){
        file.path(fish_dir, nes_path)
    })

paths_tibble <- tibble(sobj = sobj_paths, nes = nes_paths) %>%
    mutate(dataset = c("merFISH", "seqFISH")) %>%
    select(dataset, everything()) %>%
    mutate(liana = str_glue("{fish_dir}/{dataset}_liana_res.RDS"))
paths_tibble

# II) RUN LIANA and save output ----
pmap(list(paths_tibble$sobj, paths_tibble$liana), function(sobj, liana){
    message(str_glue("Loading {sobj}"))
    seurat_object <- readRDS(sobj)
    message(levels(Idents(seurat_object)))

    message("Running LIANA /w murine symbols")
    liana_res <- liana_wrap(seurat_object,
                            cellchat.params=list(organism="mouse"),
                            resource = "custom",
                            external_resource = murine_resource,
                            assay="SCT",
                            squidpy.param = list(cluster_key = "seurat_clusters"),
                            expr_prop = 0.1)
    message(levels(Idents(seurat_object)))

    message(str_glue("Saving {liana}"))
    saveRDS(liana_res, liana)

    return()
    })


# III) Fischer's Exact Test on Colocalized
fets <- pmap(.l=(list(paths_tibble$liana,
                     paths_tibble$nes,
                     paths_tibble$dataset)),
             .f=function(liana_path, nes_path, dataset_name){
                 # NES
                 nes <- readRDS(nes_path) %>%
                     as.matrix() %>%
                     reshape_coloc_estimate() %>%
                     mutate(across(c(celltype1, celltype2), ~str_replace_all(.x," ", "\\."))) %>%
                     mutate(across(c(celltype1, celltype2), ~str_replace_all(.x, "[/]", "\\.")))

                 # Load and Format LIANA res
                 liana_res <- readRDS(liana_path)
                 liana_format <- liana_res %>%
                     map(function(res) res %>%
                             mutate(across(c(ligand, receptor), str_to_title))) %>%
                     liana_aggregate() %>%
                     liana_agg_to_long()

                 liana_loc <- liana_format %>%
                     left_join(nes, by = c("source"="celltype1",
                                           "target"="celltype2"))  %>%
                     # FILTER AUTOCRINE
                     filter(source!=target) %>%
                     dplyr::mutate(localisation = case_when(estimate >= arb_thresh ~ "colocalized",
                                                            estimate < arb_thresh ~ "not_colocalized"
                     )) %>%
                     ungroup()

                 fets <- n_ranks %>%
                     map(function(n_rank){
                         run_coloc_fet(liana_loc = liana_loc,
                                       n_rank = n_rank) %>%
                             mutate(n_rank = n_rank)
                     }) %>%
                     bind_rows() %>%
                     mutate(dataset = dataset_name)
             }) %>% bind_rows()

# PLOT
boxplot_data <- fets %>%
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
    geom_hline(yintercept = 1, colour = "lightblue",
               linetype = 2, size = 0.9) +
    geom_hline(yintercept = -1, colour = "pink",
               linetype = 2, size = 0.9) +
    theme(strip.text.x = element_text(angle = 90),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
    ) +
    labs(shape=guide_legend(title="Dataset")) +
    ylab("Enrichment") +
    xlab("Number of Ranks Considered")




## FET TEST ----
nes <- readRDS(paths_tibble$nes[[2]]) %>%
    as.matrix() %>%
    reshape_coloc_estimate() %>%
    mutate(across(c(celltype1, celltype2), ~str_replace_all(.x," ", "\\."))) %>%
    mutate(across(c(celltype1, celltype2), ~str_replace_all(.x, "[/]", "\\.")))

# LIANA res
liana_res <- readRDS(paths_tibble$liana[[2]])
liana_format <- liana_res %>%
    map(function(res) res %>%
            mutate(across(c(ligand, receptor), str_to_title))) %>%
    liana_aggregate() %>%
    liana_agg_to_long()

liana_loc <- liana_format %>%
    left_join(nes, by = c("source"="celltype1",
                          "target"="celltype2"))  %>%
    # FILTER AUTOCRINE
    filter(source!=target) %>%
    dplyr::mutate(localisation = case_when(estimate >= arb_thresh ~ "colocalized",
                                           estimate < arb_thresh ~ "not_colocalized"
    )) %>%
    ungroup()

# FET
fets <- n_ranks %>%
    map(function(n_rank){
        run_coloc_fet(liana_loc = liana_loc,
                      n_rank = n_rank) %>%
            mutate(n_rank = n_rank)
    }) %>%
    bind_rows() %>% arrange(desc(enrichment))

paths_tibble

# feature space matters!!!
# SlideSeq -> show how using lower feature space destroys the results


############################# Breast Cancer ####################################
arb_thresh <- 1.645
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
tibble_dict %<>% bind_rows(tibble_dict %>% mutate(cluster_key = "celltype_minor"))


# Get Deconvolution results
deconv_results <- tibble_dict %>%
    mutate(deconv_results = pmap(tibble_dict, function(slide_name, slide_subtype, cluster_key){
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
            reshape_coloc_estimate(z_scale = TRUE) %>%
            mutate(estimate = replace_na(estimate, 0))

        return(deconv_corr_long)
    }))


# load LIANA results and get coloc
coloc = pmap(.l = list(deconv_results$slide_name,
                     deconv_results$slide_subtype,
                     deconv_results$cluster_key,
                     deconv_results$deconv_results # deconvolution results
                     ),
           .f = function(slide_name,
                         slide_subtype,
                         cluster_key,
                         deconv){
               # deconvolution subdirectory
               deconv_directory <-
                   file.path(brca_dir,
                             "deconv",
                             str_glue("{slide_subtype}_{cluster_key}"))

               # load liana
               liana_res <-
                   readRDS(file.path(deconv_directory,
                                     str_glue("{slide_subtype}_{cluster_key}_liana_res.RDS")))

               # Format LIANA
               liana_format <- liana_res %>%
                   liana_aggregate() %>%
                   liana_agg_to_long()

               liana_loc <- liana_format %>%
                   left_join(deconv, by = c("source"="celltype1",
                                            "target"="celltype2")) %>%
                   # FILTER AUTOCRINE
                   filter(source!=target) %>%
                   dplyr::mutate(localisation = case_when(estimate >= arb_thresh ~ "colocalized",
                                                          estimate < arb_thresh ~ "not_colocalized"
                                                          )) %>% ungroup()
           })

fet_tibble <- deconv_results %>%
    mutate(liana_loc = coloc)
saveRDS(fet_tibble, "data/output/spatial_out/Wu_etal_2021_BRCA/fet_tibble.RDS")

## Fishers Exact Test
fet_tibble <- readRDS("data/output/spatial_out/Wu_etal_2021_BRCA/fet_tibble.RDS")

# Define ranks and celltype for FET
cluster_key = "celltype_major"
n_ranks =
    c(#50, 100,
        # 250,
        500,
        1000, 5000,
        10000, 50000
    )
fet_tibble %<>%
    filter(cluster_key == !!cluster_key)


fet_tibble <- fet_tibble %>%
    mutate(fet_res = liana_loc %>%
               map(function(loc){
                   n_ranks %>%
                       map(function(n_rank){
                           run_coloc_fet(liana_loc = loc,
                                         n_rank = n_rank) %>%
                               mutate(n_rank = n_rank)
                       })
               }))
#
# fet_tibble <- fet_tibble %>%
#     select(-c(deconv_results,liana_loc)) %>%
#     unnest(fet_res) %>%
#     unnest(fet_res) %>%
#     unite(slide_subtype, slide_name, col="dataset")

boxplot_data <- fet_tibble %>%
    select(-c(deconv_results, liana_loc)) %>%
    unnest(fet_res) %>%
    unnest(fet_res) %>%
    unite(slide_subtype, slide_name, col="dataset") %>%
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
    geom_hline(yintercept = 1, colour = "lightblue",
               linetype = 2, size = 0.9) +
    geom_hline(yintercept = -1, colour = "pink",
               linetype = 2, size = 0.9) +
    theme(strip.text.x = element_text(angle = 90),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
    ) +
    labs(shape=guide_legend(title="Dataset")) +
    ylab("Enrichment") +
    xlab("Number of Ranks Considered")


