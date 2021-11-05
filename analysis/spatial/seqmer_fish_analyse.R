# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)

source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")


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


# III) Get LR LIANA-Colocalized
lr_coloc <- pmap(.l=(list(paths_tibble$liana,
                          paths_tibble$nes,
                          paths_tibble$dataset)
),
.f=function(liana_path, nes_path, dataset_name){
    # NES
    nes <- readRDS(nes_path) %>%
        as.matrix() %>%
        reshape_coloc_estimate() %>% # already z-scaled
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
        ungroup() %>%
        mutate(dataset = dataset_name)
}) %>%
    bind_rows()
# save liana LR score-colocalizations
saveRDS(lr_coloc, "data/output/spatial_out/seqfish_lrcoloc.RDS")
