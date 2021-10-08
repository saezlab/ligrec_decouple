library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
require(liana)
require(tidyverse)
require(magrittr)

## 1. Prerequisites
# load sc
path_to_data <- system.file(package = "SPOTlight")
cortex_sc <- readRDS(glue::glue("{path_to_data}/allen_cortex_dwn.rds"))

# load spatial
if (! "stxBrain" %in% SeuratData::AvailableData()[, "Dataset"]) {
    # If dataset not downloaded proceed to download it
    SeuratData::InstallData("stxBrain")
}
# 10X Genomics Visium Mouse Brain Dataset
anterior <- SeuratData::LoadData("stxBrain", type = "anterior1")
Seurat::SpatialDimPlot(anterior)


# prep sc
set.seed(123)
cortex_sc <- Seurat::SCTransform(cortex_sc, verbose = FALSE) %>%
    Seurat::RunPCA(., verbose = FALSE) %>%
    Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)
Seurat::DimPlot(cortex_sc,
                group.by = "subclass",
                label = TRUE) + Seurat::NoLegend()

# Get Markers
Seurat::Idents(object = cortex_sc) <- cortex_sc@meta.data$subclass
cluster_markers_all <- Seurat::FindAllMarkers(object = cortex_sc,
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = cluster_markers_all,
        file = "data/spotligh_test.RDS")
cluster_markers_all <- readRDS("data/spotligh_test.RDS")


# Run Spotlight ----
set.seed(123)
spotlight_ls <- spotlight_deconvolution(
    se_sc = cortex_sc,
    counts_spatial = anterior@assays$Spatial@counts,
    clust_vr = "subclass", # Variable in sc_seu containing the cell-type annotation
    cluster_markers = cluster_markers_all, # Dataframe with the marker genes
    cl_n = 100, # number of cells per cell type to use
    hvg = 3000, # Number of HVG to use
    ntop = NULL, # How many of the marker genes to use (by default all)
    transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
    method = "nsNMF", # Factorization method
    min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
    )
saveRDS(object = spotlight_ls, file = here::here("data/spotlight_ls.rds"))

# load spotlight deconvolution check
spotlight_ls <- readRDS("data/spotlight_ls.rds")
cluster_markers_all <- readRDS("data/spotligh_test.RDS")

nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

# Before even looking at the decomposed spots we can gain insight on how well the model performed by looking at the topic profiles for the cell types.
# The first thing we can do is look at how specific the topic profiles are for each cell type.
h <- NMF::coef(nmf_mod[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(
    h = h,
    train_cell_clust = nmf_mod[[2]])

# Check topics
topic_profile_plts[[2]]


# Remove cell types not predicted to be on the tissue
cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
decon_mtrx_sub <- decon_mtrx[, cell_types_all]
decon_mtrx_sub <- decon_mtrx_sub[, colSums(decon_mtrx_sub) > 0]

# Compute correlation
decon_cor <- cor(decon_mtrx_sub)

# Compute correlation P-value
p.mat <- corrplot::cor.mtest(mat = decon_mtrx_sub, conf.level = 0.95)

# Visualize
ggcorrplot::ggcorrplot(
    corr = decon_cor,
    p.mat = p.mat[[1]],
    hc.order = TRUE,
    type = "full",
    insig = "blank",
    lab = TRUE,
    outline.col = "lightgrey",
    method = "square",
    # colors = c("#4477AA", "white", "#BB4444"))
    colors = c("#6D9EC1", "white", "#E46726"),
    title = "Predicted cell-cell proportion correlation",
    legend.title = "Correlation\n(Pearson)") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(size = 22, hjust = 0.5, face = "bold"),
        legend.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 15),
        axis.text.x = ggplot2::element_text(angle = 90),
        axis.text = ggplot2::element_text(size = 18, vjust = 0.5))


### LIANA on cortex -----
require(liana)
require(tidyverse)
require(magrittr)
require("biomaRt")

#' Basic function to convert human to mouse genesymbols (temporary solution)
#' @param op_resource omnipath_resource as obtained via `liana::select_resource`
#'
#' @details adapted from https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
convert_to_murine <- function(op_resource){

    # query biomaRt databases
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    # obtain tibble with human and murine genesymbol
    symbols_tibble <- getLDS(attributes = c("hgnc_symbol"),
                             filters = "hgnc_symbol",
                             values = union(op_resource$source_genesymbol,
                                            op_resource$target_genesymbol),
                             mart = human,
                             martL = mouse,
                             attributesL = c("mgi_symbol")) %>%
        dplyr::rename(human_symbol = HGNC.symbol,
                      murine_symbol = MGI.symbol) %>%
        as_tibble()

    # intentionally we introduce duplicates, if needed
    # these should be resolved when LIANA is called
    # as inappropriately matched genes will not be assigned any values
    op_resource %>%
        left_join(symbols_tibble, by=c("target_genesymbol"="human_symbol")) %>%
        mutate(target_genesymbol = murine_symbol, .keep = "unused") %>%
        left_join(symbols_tibble, by=c("source_genesymbol"="human_symbol")) %>%
        mutate(source_genesymbol = murine_symbol, .keep = "unused") %>%
        filter(!is.na(target_genesymbol) | !is.na(source_genesymbol)) %>%
        filter(!is.na(target_genesymbol)) %>%
        filter(!is.na(source_genesymbol))
}

# Run LIANA ----
# Format for LIANA
cortex_sc@meta.data %<>%
    mutate(subclass = str_replace_all(subclass, "[/]", ".")) %>%
    mutate(subclass = as.factor(subclass))
Idents(cortex_sc) <- cortex_sc@meta.data$subclass

cortex_sc <- subset(cortex_sc, cells = rownames(cortex_sc@meta.data))
Seurat::DefaultAssay(cortex_sc) <- "RNA"
cortex_sc %<>% NormalizeData()
GetAssayData(cortex_sc)


# Run
op_resource <- select_resource("OmniPath")[[1]] %>%
    convert_to_murine()

liana_res <- liana_wrap(cortex_sc,
                        cellchat.params=list(organism="mouse"),
                        resource = "custom",
                        external_resource = op_resource,
                        expr_prop=0.1)

# squidpy sets gene names to upper (in the processing), revert this to title (i.e. murine)
liana_res$squidpy %<>%
    mutate_at(.vars = c("ligand", "receptor"), str_to_title)
saveRDS(liana_res, "data/spotlight_liana.rds")



# Spatial-Corr ROC here ----
# load deconv results
spotlight_ls <- readRDS("data/spotlight_ls.rds")

# corr on deconv
decon_mtrx <- spotlight_ls[[2]]
decon_cor <- cor(decon_mtrx)

# Load Liana
liana_res <- readRDS("data/spotlight_liana.rds")
liana_res %<>% liana_aggregate()


corr_response <- decon_cor %>%
    reshape2::melt() %>%
    as_tibble() %>%
    dplyr::rename(celltype1 = Var1,
                  celltype2 = Var2,
                  estimate = value) %>%
    dplyr::mutate(response = case_when(estimate >= 0.2 ~ 1,
                                       estimate < 0.2 ~ 0)) %>%
    dplyr::mutate(response = if_else(celltype1==celltype2, 0, response))

# Prep for Roc
liana_response <- liana_res %>%
    mutate(across(.cols = c(source, target), ~gsub(" ", "\\.", .x))) %>%
    filter(source!=target) %>% # REMOVE autocrine
    left_join(corr_response, by = c("source"="celltype1",
                                    "target"="celltype2")) %>%
    unite(source, ligand, target, receptor, col = "interaction") %>%
    dplyr::select(interaction, ends_with("rank"), response) %>%
    pivot_longer(cols = -c(interaction, response),
                 names_to = "method_name",
                 values_to = "predictor")%>%
    mutate(response = factor(response, levels = c(1, 0))) %>%
    mutate(predictor = predictor*-1)


# Roc
lr_space_res <- liana_response %>%
    group_by(method_name) %>%
    filter(!(method_name %in% c("mean_rank", "median_rank"))) %>%
    group_nest(.key = "corr_rank") %>%
    mutate(roc = map(corr_rank, function(df) calc_curve(df, curve="ROC"))) %>%
    mutate(prc = map(corr_rank, function(df) calc_curve(df, curve = "PR",
                                                        downsampling = TRUE,
                                                        times = 100)))
# lr_space_auroc <-
lr_space_res %>%
    dplyr::select(method_name, roc) %>%
    unnest(roc) %>%
    dplyr::select(method_name, auc) %>%
    distinct()

ggplot(lr_space_res %>%
           dplyr::select(method_name, roc) %>%
           unnest(roc), aes(x = 1-specificity,
                            y = sensitivity,
                            colour = .data$method_name)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")

# lr_space_prroc <-
lr_space_res %>%
    dplyr::select(method_name, prc) %>%
    unnest(prc) %>%
    dplyr::select(method_name, auc) %>%
    distinct()

# Correlation between score ranks and cell-cell proportion correlations -----
# cell-cell interaction proportions
n_rank <- 500
liana_res <- readRDS("data/spotlight_liana.rds")
liana_res %<>% liana_aggregate(cap = n_rank)

# load spotlight deconvolution check
spotlight_ls <- readRDS("data/spotlight_ls.rds")
decon_mtrx <- spotlight_ls[[2]]

liana_prop <- liana_res %>%
    select(source, target, ends_with("rank"), -c(mean_rank, median_rank)) %>%
    mutate(aggregate_rank = min_rank(aggregate_rank)) %>%
    pivot_longer(-c(source,target),
                 names_to = "method_name",
                 values_to = "rank") %>%
    filter(rank < n_rank) %>%
    group_by(method_name) %>%
    mutate(tot_n = n()) %>%
    group_by(source, target, method_name) %>%
    mutate(n = n()) %>%
    group_by(source, target, method_name) %>%
    mutate(prop = n/tot_n) %>%
    select(source, target, method_name, prop) %>%
    distinct() %>%
    group_by(method_name) %>%
    mutate(check_prop = sum(prop))

corr_deconv_mat <- cor(decon_mtrx) %>%
    reshape2::melt() %>%
    as_tibble() %>%
    dplyr::rename(celltype1 = Var1,
                  celltype2 = Var2,
                  estimate = value)

liana_corr <- liana_prop %>%
    left_join(corr_deconv_mat, by = c("source"="celltype1",
                                    "target"="celltype2")) %>%
    filter(source!=target) %>%
    group_by(method_name) %>%
    group_nest() %>%
    mutate(correlation = data %>%
               map(function(d) d %>%
                       mutate(estimate = replace_na(estimate, 0)) %>%
                       summarise(correlation = cor(prop,
                                                   estimate,
                                                   method = "pearson")))) %>%
    unnest(correlation)
liana_corr



#' @title Calculate AUROC and PRROC from rank-adt tibbles- `prepare_for_roc`-formated
#' `adt_rank` elements
#'
#' @param df run_benchmark roc column provided as input
#' @param downsampling logical flag indicating if the number of Negatives
#'    should be downsampled to the number of Positives
#' @param times integer showing the number of downsampling
#' @param curve whether to return a Precision-Recall Curve ("PR") or ROC ("ROC")
#' @param seed An integer to set the RNG state for random number generation. Use
#'    NULL for random number generation.
#'
#' @return tidy data frame with precision, recall, auc, n, cp, cn and coverage
#'    in the case of PR curve; or sensitivity and specificity, auc, n, cp, cn
#'    and coverage in the case of ROC.
#' @import yardstick
#'
#' @export
calc_curve = function(df,
                      downsampling = FALSE,
                      times = 1000,
                      curve = "ROC",
                      seed = 420){
    set.seed(seed)

    if(curve=="PR"){
        res_col_1 <- "precision"
        res_col_2 <- "recall"
        curve_fun = yardstick::pr_curve
        auc_fun = yardstick::pr_auc
    }
    else if(curve=="ROC"){
        res_col_1 <- "sensitivity"
        res_col_2 <- "specificity"
        curve_fun = yardstick::roc_curve
        auc_fun = yardstick::roc_auc
    }

    if (sum(which(df$response == 0)) == nrow(df)){
        return(as_tibble(NULL))
    }

    cn = df %>% filter(.data$response == 0)
    cp = df %>% filter(.data$response == 1)

    feature_coverage = length(unique(df$interaction))

    if(downsampling == TRUE){
        num_tp = nrow(cp)

        res = map_df(seq(from=1, to=times, by=1), function(i) {
            df_sub = sample_n(cn, num_tp, replace=TRUE) %>%
                bind_rows(cp)

            r_sub = df_sub %>%
                curve_fun(.data$response, .data$predictor)

            auc = df_sub %>%
                auc_fun(.data$response, .data$predictor) %>%
                pull(.data$.estimate)

            res_sub = tibble({{ res_col_1 }} := r_sub %>% pull(res_col_1),
                             {{ res_col_2 }} := r_sub %>% pull(res_col_2),
                             th = r_sub$.threshold,
                             auc = auc,
                             n = length(which(df$response == 1)),
                             cp = nrow(cp),
                             cn = nrow(cn),
                             coverage = feature_coverage) %>%
                mutate("run" = i)

        })
        # Get Average AUC
        res <- res %>% dplyr::rename("raw_auc" = auc)
        # auc is the mean of all iterations, raw_auc is the value per iteration
        res$auc <- sum(res$raw_auc)/length(res$raw_auc)
        res$cn <- nrow(cp)

    } else {
        r = df %>%
            curve_fun(.data$response, .data$predictor)
        auc = df %>%
            auc_fun(.data$response, .data$predictor)

        res = tibble({{ res_col_1 }} := r %>% pull(res_col_1),
                     {{ res_col_2 }} := r %>% pull(res_col_2),
                     th = r$.threshold,
                     auc = auc$.estimate,
                     n = length(which(df$response == 1)),
                     cp = nrow(cp),
                     cn = nrow(cn),
                     coverage = feature_coverage) %>%
            arrange(!!res_col_1, !!res_col_2)
    }
    return(res)
}

