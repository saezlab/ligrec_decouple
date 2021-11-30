# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)

source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")
source("analysis/comparison/comparison_utils.R")

fish_dir <- "data/input/spatial/fishes"
murine_resource <- readRDS("data/input/murine_omnipath.RDS")


########################### SeqFISH, merFISH ###################################
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


# II) RUN LIANA and save output ----
pmap(list(paths_tibble$sobj,
          paths_tibble$liana),
     function(sobj, liana){

    message(str_glue("Loading {sobj}"))
    seurat_object <- readRDS(sobj)
    message(levels(Idents(seurat_object)))

    message("Running LIANA /w murine symbols")
    liana_res <- liana_wrap(seurat_object,
                            resource = "custom",
                            external_resource = murine_resource,
                            expr_prop = 0.1,
                            cellchat.params = list(nboot=1000,
                                                   expr_prop = 0,
                                                   organism="mouse"),
                            assay="SCT",
                            method = c('call_natmi',
                                       'call_connectome', 'logfc', 'cellchat',
                                       'call_sca', 'cellphonedb', "cytotalk"
                            ))
    message(levels(Idents(seurat_object)))

    message(str_glue("Saving {liana}"))
    saveRDS(liana_res, liana)

    return()
})


# Add setting
paths_tibble %<>%
    mutate(setting = "comp_n") %>%
    bind_rows(., .)
paths_tibble$setting[c(3,4)] <- "house_n"
paths_tibble


# III) Get LR LIANA-Colocalized
# Set to Independent and Max and liana_aggregate_enh
lr_coloc <- pmap(
    .l=(list(
        paths_tibble$liana,
        paths_tibble$nes,
        paths_tibble$dataset,
        paths_tibble$setting)),
    .f=function(liana_path,
                nes_path,
                dataset_name,
                setting){
    set_aggregation_settings(setting)

    # NES
    nes <- readRDS(nes_path) %>%
        as.matrix() %>%
        reshape_coloc_estimate() %>% # already z-scaled
        mutate(across(c(celltype1, celltype2), ~str_replace_all(.x," ", "\\."))) %>%
        mutate(across(c(celltype1, celltype2), ~gsub("[/]", "\\.", .x)))

    # Load and Format LIANA res
    liana_res <- readRDS(liana_path)
    liana_format <- liana_res %>%
        map(function(res) res %>%
                mutate(across(c(ligand, receptor), str_to_title))) %>%
        liana_aggregate_enh(
            filt_de_pvals = TRUE,
            de_thresh = de_thresh, # we only filter Connectome DEs
            filt_outs = TRUE,
            pval_thresh = pval_thresh,
            sca_thresh = 0,
            .score_mode = .score_specs,
            .eval = "independent"
        ) %>%
        liana_agg_to_long()

    liana_loc <- liana_format %>%
        left_join(nes, by = c("source"="celltype1",
                              "target"="celltype2"))  %>%
        # FILTER AUTOCRINE
        filter(source!=target) %>%
        ungroup() %>%
        mutate(dataset = dataset_name) %>%
        mutate(setting = setting) %>%
        dplyr::mutate(localisation =
                          case_when(estimate >= 1.645 ~ "colocalized",
                                    estimate < 1.645 ~ "not_colocalized")
                      )
    }) %>%
    bind_rows()
# save liana LR score-colocalizations
saveRDS(lr_coloc, str_glue("data/output/spatial_out/fish/fish_lrcoloc.RDS"))
