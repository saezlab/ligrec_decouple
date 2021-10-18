# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)

source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")

# merFISH genes
# seqfish_obj <- readRDS("data/input/spatial/fishes/seqFISH_seurat.RDS")

# MOUSE BRAIN ATLAS ----
brain_dir <- "data/input/spatial/brain_cortex/"
murine_resource <- readRDS("data/input/murine_omnipath.RDS")

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
            reshape_coloc_estimate(z_scale = FALSE)

        return(deconv_corr_long)
    }) %>%
    setNames(slides)

# Bind lr and coloc
lr_coloc <- deconv_results %>%
    map2(names(.), function(deconv_cor_formatted, slide_name){
        # Assign colocalisation to liana results in long according to a threshold
        liana_loc <- liana_format %>%
            left_join(deconv_cor_formatted, by = c("source"="celltype1",
                                                   "target"="celltype2")) %>%
            # FILTER AUTOCRINE
            filter(source!=target) %>%
            ungroup() %>%
            mutate(dataset = slide_name)
    }) %>%
    bind_rows()
# save liana LR score-colocalizations
saveRDS(lr_coloc, "data/output/spatial_out/brain_lrcoloc.RDS")


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




############################# Breast Cancer ####################################
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


# NOTE WE KEEP ONLY celltype minor (major are also random) :D
tibble_dict %<>% filter(cluster_key == "celltype_minor")

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
            reshape_coloc_estimate(z_scale = FALSE) %>%
            mutate(estimate = replace_na(estimate, 0))

        return(deconv_corr_long)
    }))


# load LIANA results and get coloc
lr_coloc <- pmap(.l = list(deconv_results$slide_name,
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
                                                          )) %>%
                   ungroup() %>%
                   mutate(dataset = slide_name) %>%
                   mutate(subtype = slide_subtype)
           }) %>%
    bind_rows()
saveRDS(lr_coloc, "data/output/spatial_out/brca_lrcoloc.RDS")


############################## PLOT ############################################
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)

source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")

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
    distinct()  %>%
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
brca_lr_coloc <- brca_lr_coloc %>%
    select(-localisation) %>%
    left_join(coloc)


# 1b) Initial check
brca_lr_coloc %>%
    check_coloc()
hist(brca_lr_coloc$estimate)

# 2b) FET boxplot ON CONSENSUS
brca_lr_coloc %>%
    get_fet_boxplot_data(., n_ranks = n_ranks) %>%
    get_spatial_boxplot()


# IV) Bind ALL ----
all_lr_coloc <- bind_rows(brain_lr_coloc,
                          seqfish_lr_coloc,
                          brca_lr_coloc)
saveRDS(all_lr_coloc, "data/output/spatial_out/all_lrcoloc.RDS")

# 1) Initial Check All
all_lr_coloc %>%
    check_coloc()
summary(as.factor(all_lr_coloc$localisation))

# 2) FET boxplot All
all_lr_coloc %>%
    get_fet_boxplot_data(., n_ranks = n_ranks) %>%
    get_spatial_boxplot()

# all_lr_coloc %<>%
#     dplyr::mutate(localisation = case_when(estimate >= arb_thresh ~ "colocalized",
#                                            estimate < arb_thresh ~ "not_colocalized"
#     ))


# 3) AUROC/AUPRC Curves on ALL -----
# A) Positive Correlation AUC ----
all_lr_coloc <- readRDS("data/output/spatial_out/all_lrcoloc.RDS")
all_lr_coloc #%<>%
    # group_by(method_name, dataset) %>%
    # slice_min(prop = 0.1, order_by = predictor) %>%
    # ungroup()
     # filter to proportion (e.g. top 10% of ranks)

all_lr_coloc <- all_lr_coloc %>%
    mutate(response = if_else(localisation=="colocalized",
                              1,
                              0) %>% as.factor()) %>%
    select(-c(subtype, loc_consensus))


all_lr_auc <- all_lr_coloc %>%
    group_by(method_name, dataset) %>%
    group_nest(.key = "method_dataset")  %>%
    mutate(prc = map(method_dataset,
                     function(df)
                         calc_curve(df,
                                    curve = "PR",
                                    downsampling = TRUE,
                                    source_name = "interaction",
                                    auc_only = TRUE))) %>%
    mutate(roc = map(method_dataset,
                     function(df) calc_curve(df,
                                             curve="ROC",
                                             downsampling = FALSE,
                                             times = 100,
                                             source_name = "interaction",
                                             auc_only = TRUE))) %>%
    select(-method_dataset)
saveRDS(all_lr_auc, "data/output/spatial_out/auc_positive.RDS")
all_lr_auc



# Load results
# ROC
all_lr_auc <- readRDS("data/output/spatial_out/brain_cortex/rocs_positive.RDS")
pos_roc <- get_auroc_heat(all_lr_auc, "roc") # all random
pos_roc

# PRC
pos_prc <- get_auroc_heat(all_lr_auc, "prc")
pos_prc
gc()

### NEGATIVE AUROC - maybe works better for BRCA?
