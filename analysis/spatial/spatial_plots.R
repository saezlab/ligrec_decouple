############################## PLOT ############################################
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)

source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

# z_thresh <- 1.645 # alpha = 0.05 # 0.842 # alpha = 0.2
corr_thresh <- 0.25
n_ranks = c(50, 100,
            500, 1000,
            2500, 5000,
            10000)

# Bind ALL ----
all_lr_coloc <- bind_rows(readRDS("data/output/spatial_out/brain_cortex/coloc_inpdependent_specs.RDS") %>%
                              dplyr::mutate(
                                  localisation = case_when(estimate >= corr_thresh ~ "colocalized",
                                                           estimate < corr_thresh ~ "not_colocalized")
                                  ),
                          # seqfish_lr_coloc,
                          readRDS("data/output/spatial_out/Wu_etal_2021_BRCA/coloc_independent_specs.RDS") %>%
                              dplyr::mutate(
                                  localisation = case_when(estimate >= corr_thresh ~ "colocalized",
                                                           estimate < corr_thresh ~ "not_colocalized")
                              )
                          )
saveRDS(all_lr_coloc, "data/output/spatial_out/all_lrcoloc.RDS")

# 1) Initial Check All
all_lr_coloc <- readRDS("data/output/spatial_out/all_lrcoloc.RDS")
all_lr_coloc %>% check_coloc()
summary(as.factor(all_lr_coloc$localisation))

# 2.1) FET boxplot All by slide
# all_lr_coloc %>%
#     get_fet_boxplot_data(., n_ranks = n_ranks) %>%
#     get_spatial_boxplot()


# 2.2) FET boxplot All by dataset type
type_lr_coloc <- all_lr_coloc %>%
    filter(!(dataset %in% c("merFISH", "seqFISH"))) %>% # Remove FISH
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
saveRDS(type_lr_coloc, "data/output/spatial_out/all_fets_specs.RDS")


# Boxplot
cairo_pdf(file.path("data/output/temp/",
                    str_glue("spatial_fets_specs.RDS")),
          height = 24,
          width = 32,
          family = 'DINPro')
print(
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
)
dev.off()



