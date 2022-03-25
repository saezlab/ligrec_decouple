require(tidyverse)
require(magrittr)
require(liana)
require(SingleCellExperiment)
require(scuttle)
require(decoupleR)

source("scripts/utils.R")
source("scripts/liana_cytosig.R")
ligrec_path <- "/media/dbdimitrov/SSDDimitrov/Repos/ligrec_decouple/"

source(file.path(ligrec_path, "analysis/spatial/spatial_src.R"))
source(file.path(ligrec_path, "src/eval_utils.R"))
source(file.path(ligrec_path, "src/plot_utils.R"))

# Geenrate PM-only resource to which we filter
pm_omni <- generate_omni(
    loc_consensus_percentile = 51, # increase localisation consensus threshold
    consensus_percentile = NULL,
    # include only PM-bound proteins
    transmitter_topology = c(
        'plasma_membrane_transmembrane',
        'plasma_membrane_peripheral'
    ),
    receiver_topology = c(
        'plasma_membrane_transmembrane',
        'plasma_membrane_peripheral'
    ),
    min_curation_effort = 1,
    ligrecextra = FALSE,
    remove_complexes = FALSE, # keep complexes
    simplify = TRUE # do simplify
) %>%
    unite(col = "key", remove = FALSE, source_genesymbol, target_genesymbol) %>%
    # Note ideal but a quick fix for murine
    mutate(across(c("source_genesymbol", "target_genesymbol"), ~toupper(.x)))


#
conditions <- c("ER", "TNBC", "brain")


#
fet_results <- map(conditions, function(condition){

    # load aggregated liana + spatial correlations
    liana_agg <- readRDS(file.path(ligrec_path,
                                   "data/output/aggregates/",
                                   str_glue("{condition}_max_mixed_liana_res.RDS")))

    spatial_corr_path <- file.path(ligrec_path, str_glue("data/output/spatial_out/deconv_summ/{condition}_coloc_corr.RDS"))

    liana_agg_filt <- liana_agg %>%
        # Generate a key by which we filter down to PM-bound only
        mutate(across(c("ligand", "receptor"), ~toupper(.x))) %>%
        unite(col="key", ligand, receptor, remove = FALSE) %>%
        # filter(key %in% pm_omni$key) %>% # SIMPLY CHANGE ! to be SURFACE or SECRETED
        select(source, target, ligand, receptor, ends_with("rank"), starts_with("cc")) %>%
        na.omit() %>%
        mutate(across(ends_with("rank"), ~min_rank(.x)))

    saveRDS(liana_agg_filt, str_glue("data/output/liana_filt_{condition}.RDS"))

    fet_coloc <- get_lr_colocalized(str_glue("data/output/liana_filt_{condition}.RDS"),
                                    spatial_corr_path,
                                    condition = condition,
                                    corr_thresh = 1.645,
                                    n_ranks = c(50, 100, 250,
                                                500, 1000,
                                                2500, 5000, 10000
                                    ))
    return(fet_coloc)
    }) %>%
    setNames(conditions) %>%
    bind_rows()
saveRDS(fet_results, "data/output/fet_spatial_only_all.RDS")

fet_results_pm <- readRDS("data/output/fet_spatial_only_pm.RDS")
fet_results_secreted <- readRDS("data/output/fet_spatial_only_secreted.RDS")
fet_results_all <- readRDS("data/output/fet_spatial_only_all.RDS")

ggplot(fet_results_secreted,
       aes(x = n_rank, y = odds_ratio,
           color = dataset, group=dataset)) +
    geom_point(aes(shape = dataset), size = 3, alpha = 0.8) +
    facet_grid(~method_name, scales='free_x', space='free', switch="x") +
    theme_bw(base_size = 21) +
    geom_line(size = 1.5, alpha = 0.6) +
    geom_hline(yintercept = 1, colour = "lightblue",
               linetype = 2, size = 1.7) +
    theme(strip.text.x = element_text(angle = 90, face="bold", colour="white"),
          axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5),
          strip.background = element_rect(fill="navyblue"),
          legend.title = element_text(size = 23),
          legend.text = element_text(size = 21)
    ) +
    labs(colour=guide_legend(title="Dataset")) +
    ylab("Odds Ratio") +
    xlab("Ranked Interactions Range") +
    guides(shape = "none")




