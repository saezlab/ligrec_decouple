# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)
require(SingleCellExperiment)
require(decoupleR)

source("analysis/cytosig/cytosig_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

source("analysis/comparison/comparison_utils.R")
set_aggregation_settings("comp_n")
.eval = "independent"
score_mode = "comp_n"
generate = TRUE
path_tibble = NULL
outpath = "data/output/cytosig_out/"

liana_res <- readRDS("data/output/comparison_out/crc_liana_res.RDS") %>%
    transpose() %>%
    pluck("OmniPath") %>%
    liana_aggregate_enh(sca_thresh = sca_thresh,
                        pval_thresh = pval_thresh,
                        de_thresh = de_thresh,
                        .eval = "max",
                        .score_mode =.score_specs
                        )
saveRDS(liana_res, "data/output/temp/crc_aggregated.RDS")
liana_res <- readRDS("data/output/temp/crc_aggregated.RDS")


path_tibble %<>% `%||%`(
    tibble(dataset = c("CRC"),
           # BRCA seurat objects used throughout the manuscript
           seurat_path = c(
               # "data/input/spatial/Wu_etal_2021_BRCA/deconv/ER_celltype_minor/ER_celltype_minor_seurat.RDS",
               # "data/input/spatial/Wu_etal_2021_BRCA/deconv/HER2_celltype_minor/HER2_celltype_minor_seurat.RDS",
               # "data/input/spatial/Wu_etal_2021_BRCA/deconv/TNBC_celltype_minor/TNBC_celltype_minor_seurat.RDS"
               "data/input/comparison/crc_data/crc_korean_form.rds"
           ),
           # liana results generated from extract_evals.sh
           liana_path = c(
               "data/output/temp/crc_aggregated.RDS"
               # file.path("data/output/aggregates", str_glue("ER_{.eval}_{score_mode}_liana_res.RDS")),
               # file.path("data/output/aggregates", str_glue("HER2_{.eval}_{score_mode}_liana_res.RDS")),
               # file.path("data/output/aggregates", str_glue("TNBC_{.eval}_{score_mode}_liana_res.RDS"))
           ))

)
path_tibble

# get cytosig
cytosig_net <- load_cytosig()


seurat_object <- readRDS("data/input/comparison/crc_data/crc_korean_form.rds")
pseudo <- get_pseudobulk(seurat_object,
                         assay = "RNA",
                         expr_prop = 0.1,
                         sum_count_thresh = 5)
message("Pseudobulk Cytokine Enrichment")


# Cytokine Activity Enrichment
pseudo_cytosig <- pseudo %>%
    mutate(cytosig_res = logcounts %>%
               map(function(logc){
                   run_mlm(
                       logc,
                       cytosig_net,
                       .source = "cytokine",
                       .target = "target",
                       .mor = "mor",
                       .likelihood = "weight",
                       sparse = FALSE
                   ) %>%
                       # rename
                       select(cytokine=source,
                              NES=score,
                              p_value)
               })
    )  %>%
    select(celltype, cytosig_res)



### Join LIANA to cytosig
# Cytosig: add aliases (as in OP) and cytokine family members if appropriate
aliases_and_families <- list("CD40L" = "CD40LG",
                             "GSFC" = "CSF3",
                             "IFN1" = c("IFNA1", "IFNA2", "IFNA10",
                                        "IFNA7", "IFNA21", "IFNA5",
                                        "IFNA14", "IFNA17", "IFNA6",
                                        "IFNA4",  "IFNA16", "IFNA8" ),
                             "IFNL" = c("IFNL1", "IFNL2", "IFNL3", "IFNL4"),
                             "IL12" = c("IL12A", "IL12B"),
                             "IL36" = c("IL36A", "IL36B", "IL36G", "IL36RN"),
                             "MCSF" = "CSF1",
                             "TNFA" = c("TNF", "TNFA"),
                             "TRAIL" = c("TNFSF10"),
                             "TWEAK" = "TNFSF12")
alias_tib <- tibble(cytokine = names(aliases_and_families),
                    aliases = aliases_and_families) %>%
    unnest(aliases)

# cytosig results format
cytosig_res <- pseudo_cytosig %>%
    unnest(cytosig_res) %>%
    left_join(alias_tib, by = "cytokine") %>%
    mutate(aliases = if_else(is.na(aliases),
                             cytokine,
                             aliases)) %>%
    select(-cytokine, cytokine=aliases) %>%
    select(cytokine, celltype, NES, p_value) %>%
    { if(z_scale) z_transform_nes(.)  else . } # cluster-specific or not


# Format and join cytokine activities and prep for AUC
liana_cytosig <- liana_res %>%
    # keep only interactions in which the ligand is a cytokine present in cytosig
    filter(ligand %in% unique(cytosig_res$cytokine)) %>%
    dplyr::select(ligand, target, ends_with("rank"), -c(mean_rank, median_rank)) %>%
    mutate(aggregate_rank = min_rank(aggregate_rank)) %>%
    # here we join by ligand from liana and cytokine in cytosig,
    # as well as target cluster from liana and the celltype for which we predicted
    # the ligand activities
    left_join(cytosig_res, by=c("ligand"="cytokine", "target"="celltype")) %>%
    pivot_longer(ends_with("rank"), names_to = "method_name", values_to = "predictor") %>%
    mutate(predictor = predictor*-1) %>% # convert rankings
    unite(ligand, target, col = "cytokine_in_target")


# Prep for ROC
cytosig_eval <- liana_cytosig %>%
    # correct p
    mutate(adj_pvalue = p.adjust(p_value, "fdr")) %>%
    mutate(response = if_else(NES > NES_thresh & adj_pvalue <=0.05,
                              1,
                              0)) %>%
    mutate(response = factor(response, levels = c(1, 0))) # first level is the truth
print(cytosig_eval %>% arrange(NES))

message("Calculating AUCs")
# ROC and Correlations
cytosig_eval %<>%
    group_by(method_name) %>%
    group_nest(.key = "cyto_liana") %>%
    mutate(roc = map(cyto_liana,
                     function(df) calc_curve(df,
                                             curve="ROC",
                                             downsampling = FALSE,
                                             source_name = "cytokine_in_target",
                                             auc_only = TRUE)
    )) %>%
    mutate(prc = map(cyto_liana,
                     function(df) calc_curve(df,
                                             curve="PR",
                                             downsampling = TRUE,
                                             times = 100,
                                             source_name = "cytokine_in_target",
                                             auc_only = TRUE))) %>%
    mutate(corr = cyto_liana %>%
               map(function(df){
                   cor.test(df[["NES"]],
                            df[["predictor"]],
                            method="spearman",
                            exact = FALSE) %>%
                       broom::tidy()
               }))



# Run Cytosig Evaluation
cytosig_eval <- path_tibble %>%
    mutate(cytosig_res =
               pmap(path_tibble,
                    function(dataset, seurat_path, liana_path){
                        # Load Seurat Object
                        message(str_glue("Now running: {dataset}"))

                        seurat_object <- readRDS(seurat_path)
                        print(seurat_object)

                        # read liana
                        liana_res <- readRDS(liana_path)

                        # Run Cytosig Evaluation
                        cyto_res <- run_cytosig_eval(seurat_object = seurat_object,
                                                     liana_res = liana_res,
                                                     cytosig_net = cytosig_net,
                                                     z_scale = FALSE,
                                                     expr_prop = 0.1,
                                                     assay = "RNA",
                                                     sum_count_thresh = 5,
                                                     NES_thresh = 0,
                                                     subtype = dataset,
                                                     generate = generate) #!
                        gc()
                        return(cyto_res)
                    }))



saveRDS(cytosig_eval,
        file.path(outpath,
                  str_glue("cytosig_res_{.eval}_{score_mode}.RDS")))
