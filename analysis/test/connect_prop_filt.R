require(liana)
require(tidyverse)
require(magrittr)
require(Seurat)
require(SingleCellExperiment)
source("analysis/citeseq/citeseq_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

# Input
# liana_path <- system.file(package = "liana")
# seurat_object <-
#     readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
seurat_object <- readRDS("data/input/citeseq/5k_pbmcs_nextgem/5k_pbmcs_nextgem_seurat.RDS")
seurat_object@meta.data %<>% mutate(seurat_clusters = as.factor(seurat_clusters))
Idents(seurat_object) <- seurat_object@meta.data$seurat_clusters

# Run liana
liana_res <- liana_wrap(seurat_object,
                        resource = c('OmniPath'),
                        method = c("natmi", "connectome", "logfc",
                                   "sca", "scconnect", "cytotalk"
                                   #, "squidpy", "call_natmi", "cellphonedb",  "cellchat"
                                   ),
                        expr_prop = 0.1)
saveRDS(liana_res, "data/output/temp/small_test2.RDS")

# conn results (UNFILT)
liana_res <- readRDS("data/output/temp/small_test2.RDS")

liana_res_full <- liana_res %>% liana_aggregate()
liana_res_intersect <- liana_res_full %>%
    mutate(across(ends_with("rank"), ~ifelse(.x==max(.x), NA, .x))) %>%
    na.omit()
liana_res_ind <- liana_res_full %>%
    mutate(across(ends_with("rank"), ~ifelse(.x==max(.x), NA, .x)))


# Citeseq check ----
# Get Symbols of Receptors from OP
op_resource <- select_resource("OmniPath")[[1]]
murine_resource <- readRDS("data/input/murine_omnipath.RDS")
arbitrary_thresh = 1.645 # one-tailed alpha = 0.05
organism <- "human"
receptor_syms <- str_to_symbol(op_resource$target_genesymbol, organism)

# Get ADT Symbols
adt_symbols <- rownames(seurat_object@assays$ADT)

# Obtain adt genesymbols
# obtain both mouse and human but convert to title
adt_aliases <- get_adt_aliases(adt_symbols = adt_symbols,
                               organism = "mouse") %>%
    bind_rows(get_adt_aliases(adt_symbols = adt_symbols,
                              organism = "human")) %>%
    mutate(across(everything(), str_to_title)) %>%
    distinct()
adt_aliases <- get_adt_aliases(adt_symbols = adt_symbols,
                               organism = organism)


# ROC FULL
adt_lr_roc <-
    generate_specificity_roc(
        seurat_object = seurat_object,
        liana_res = liana_res_full,
        op_resource = op_resource,
        cluster_key = "seurat_clusters",
        organism=organism,
        receptor_syms = receptor_syms,
        adt_symbols = adt_symbols,
        adt_aliases = adt_aliases,
        arbitrary_thresh = arbitrary_thresh)

adt_lr_roc %>%
    dplyr::select(method_name, roc) %>%
    unnest(roc) %>%
    dplyr::select(method_name, auc) %>%
    distinct()

# Intersect
adt_lr_roc_intersect <-
    generate_specificity_roc(
        seurat_object = seurat_object,
        liana_res = liana_res_intersect,
        op_resource = op_resource,
        cluster_key = "seurat_clusters",
        organism=organism,
        receptor_syms = receptor_syms,
        adt_symbols = adt_symbols,
        adt_aliases = adt_aliases,
        arbitrary_thresh = arbitrary_thresh)

# Independent
adt_lr_roc_ind <-
    generate_specificity_roc(
        seurat_object = seurat_object,
        liana_res = liana_res_ind,
        op_resource = op_resource,
        cluster_key = "seurat_clusters",
        organism=organism,
        receptor_syms = receptor_syms,
        adt_symbols = adt_symbols,
        adt_aliases = adt_aliases,
        arbitrary_thresh = arbitrary_thresh)

# Full
adt_lr_roc %>%
    dplyr::select(method_name, prc) %>%
    unnest(prc) %>%
    dplyr::select(method_name, auc) %>%
    distinct()

# Intersect
adt_lr_roc_intersect %>%
    dplyr::select(method_name, prc) %>%
    unnest(prc) %>%
    dplyr::select(method_name, auc) %>%
    distinct()

# Partial
adt_lr_roc_part %>%
    dplyr::select(method_name, prc) %>%
    unnest(prc) %>%
    dplyr::select(method_name, auc) %>%
    distinct()


# Independent
adt_lr_roc_ind %>%
    dplyr::select(method_name, prc) %>%
    unnest(prc) %>%
    dplyr::select(method_name, auc) %>%
    distinct()



