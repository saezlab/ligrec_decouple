############################## PLOT ############################################
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)

source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

z_thresh <- 1.645 # alpha = 0.05 # 0.842 # alpha = 0.2
corr_thresh <- 0.25
n_ranks = c(50, 100,
            500, 1000,
            2500, 5000,
            10000)


### I) Mouse BRAIN visium ----
# 1) Initiall Check
brain_lr_coloc <- readRDS("data/output/spatial_out/brain_lrcoloc.RDS")

# assign coloc threshold
brain_lr_coloc %<>%
    dplyr::mutate(localisation = case_when(estimate >= corr_thresh ~ "colocalized",
                                           estimate < corr_thresh ~ "not_colocalized"))

hist(brain_lr_coloc$estimate)
# Check
brain_lr_coloc %>%
    check_coloc()

# 2) FET boxplot
brain_lr_coloc %>%
    get_fet_boxplot_data(., n_ranks = n_ranks) %>%
    get_spatial_boxplot()



### II) Mouse BRAIN Seq/MerFish -----------------------------------------------
### Seq/MerFISH
seqfish_lr_coloc <- readRDS("data/output/spatial_out/seqfish_lrcoloc.RDS") %>%
    dplyr::mutate(localisation = case_when(estimate >= z_thresh ~ "colocalized",
                                           estimate < z_thresh ~ "not_colocalized"))


# 1) Initial Check
hist(seqfish_lr_coloc$estimate)
seqfish_lr_coloc %>%
    check_coloc()

# 2) FET boxplot
seqfish_lr_coloc %>%
    get_fet_boxplot_data(., n_ranks = n_ranks) %>%
    get_spatial_boxplot()


### III) BRCA ------
brca_lr_coloc <- readRDS("data/output/spatial_out/brca_lrcoloc.RDS")

# A) Coloc by individual slides
# establish colocalisation truth
brca_lr_coloc %<>%
    dplyr::mutate(localisation = case_when(estimate >= corr_thresh ~ "colocalized",
                                           estimate < corr_thresh ~ "not_colocalized"))


# 1a) Initial check
brca_lr_coloc %>%
    check_coloc()
hist(brca_lr_coloc$estimate)

# 2a) FET boxplot ON CONSENSUS
brca_lr_coloc %>%
    get_fet_boxplot_data(., n_ranks = n_ranks) %>%
    get_spatial_boxplot()


# B) Consensus by Cancer Subtype
coloc <- brca_lr_coloc %>%
    select(source, target, dataset, subtype, localisation) %>%
    distinct() %>%
    mutate(loc_consensus = if_else(localisation=="colocalized",
                                   1,
                                   0)) %>%
    group_by(source, target, subtype)  %>%
    mutate(loc_consensus = sum(loc_consensus)) %>%
    distinct() %>%
    arrange(desc(loc_consensus)) %>%
    mutate(localisation = if_else(subtype == "TNBC" & loc_consensus>=3,
                                  "colocalized",
                                  "not_colocalized")) %>%
    mutate(localisation = if_else(subtype == "ER" & loc_consensus==2,
                                  "colocalized",
                                  localisation
    ))
brca_lr_cons <- brca_lr_coloc %>%
    select(-localisation) %>%
    left_join(coloc)


# 1b) Initial check
brca_lr_cons %>%
    check_coloc()
hist(brca_lr_cons$estimate)

# 2b) FET boxplot ON CONSENSUS
brca_lr_cons %>%
    get_fet_boxplot_data(., n_ranks = n_ranks) %>%
    get_spatial_boxplot()


# IV) Bind ALL ----
all_lr_coloc <- bind_rows(brain_lr_coloc,
                          # seqfish_lr_coloc,
                          brca_lr_coloc)
saveRDS(all_lr_coloc, "data/output/spatial_out/all_lrcoloc.RDS")

# 1) Initial Check All
all_lr_coloc <- readRDS("data/output/spatial_out/all_lrcoloc.RDS")
all_lr_coloc %>% check_coloc()
summary(as.factor(all_lr_coloc$localisation))

# 2.1) FET boxplot All by slide
all_lr_coloc %>%
    get_fet_boxplot_data(., n_ranks = n_ranks) %>%
    get_spatial_boxplot()


# 2.2) FET boxplot All by dataset type
type_lr_coloc <- all_lr_coloc %>%
    filter(!(dataset %in% c("merFISH", "seqFISH"))) %>% # Remove FISHes
    get_fet_boxplot_data(., n_ranks = n_ranks) %>%
    # recode slide by dataset type
    mutate(dataset_type = dplyr::recode(dataset,
                                        "Cortex Anterior 1" = "Brain Cortex",
                                        "Cortex Anterior 2" = "Brain Cortex",
                                        "Cortex Posterior 1" = "Brain Cortex",
                                        "Cortex Posterior 2" = "Brain Cortex",
                                        "TNBC1 (1142243F)" = "TNBC",
                                        "TNBC2 (1160920F)" = "TNBC",
                                        "TNBC3 (CID4465)" = "TNBC",
                                        "TNBC4 (CID44971)" = "TNBC",
                                        "ER1 (CID4290)" = "ER-positive BC",
                                        "ER2 (CID4535)" = "ER-positive BC")) %>%
    mutate(dataset = recode_datasets(dataset))

# Boxplot
ggplot(type_lr_coloc,
       aes(x = n_rank,
           y = odds_ratio,
           color = dataset_type,
           fill = dataset_type)) +
    geom_boxplot(alpha = 0.2,
                 outlier.size = 1.5,
                 width = 0.8)  +
    # geom_point(aes(shape = dataset), size = 2) +
    scale_shape_manual(values = rep(1:12, len =  length(unique(type_lr_coloc$dataset)))) +
    facet_grid(~method_name, scales='free_x', space='free', switch="x") +
    theme_bw(base_size = 24) +
    geom_hline(yintercept = 1, colour = "black",
               linetype = 2, size = 1.3) +
    theme(strip.text.x = element_text(angle = 90),
          axis.text.x = element_text(angle = 45, hjust=1)
    ) +
    # guides(color = "none") +
    labs(shape = guide_legend(title="Visium slide"),
         color = guide_legend(title="Dataset type")) +
    guides(fill = "none") +
    ylab("Odds Ratio") +
    xlab("#Ranks Considered")


# 3) AUROC/AUPRC Curves on ALL -----
# all_lr_auc <- all_lr_coloc %>%
#     group_by(method_name, dataset) %>%
#     slice_min(prop = 0.5, order_by = predictor) %>% # only look at the peak of the iceberg - too many interactions to be able to build a ROC curve
#     ungroup() %>%
#     mutate(response = if_else(localisation=="colocalized",
#                               1,
#                               0) %>% as.factor()) %>%
#     # Get AUCs
#     group_by(method_name, dataset) %>%
#     group_nest(.key = "method_dataset")  %>%
#     # mutate(prc = map(method_dataset,
#     #                  function(df)
#     #                      calc_curve(df,
#     #                                 curve = "PR",
#     #                                 downsampling = TRUE,
#     #                                 source_name = "interaction",
#     #                                 auc_only = TRUE))) %>%
#     mutate(roc = map(method_dataset,
#                      function(df) calc_curve(df,
#                                              curve="ROC",
#                                              downsampling = FALSE,
#                                              times = 100,
#                                              source_name = "interaction",
#                                              auc_only = TRUE))) %>%
#     select(-method_dataset)
# saveRDS(all_lr_auc, "data/output/spatial_out/auc_positive.RDS")
# all_lr_auc
#
#
#
# # Load results
# # ROC
# all_lr_auc <- readRDS("data/output/spatial_out/auc_positive.RDS")
# pos_roc <- get_auroc_heat(all_lr_auc, "roc",
#                           auc_min = 0, auc_max = 1) # all random
# pos_roc
#
# # PRC
# pos_prc <- get_auroc_heat(all_lr_auc, "prc")
# pos_prc
# gc()
