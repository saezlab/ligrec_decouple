library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
library(SPOTlight)
require(liana)
require(tidyverse)
require(magrittr)






## 1. Prerequisites
# load sc
path_to_data <- system.file(package = "SPOTlight")
cortex_sc <- readRDS(glue::glue("{path_to_data}/allen_cortex_dwn.rds"))


# prep sc
set.seed(123)
cortex_sc <- Seurat::SCTransform(cortex_sc, verbose = FALSE) %>%
    Seurat::RunPCA(., verbose = FALSE) %>%
    Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)
Seurat::DimPlot(cortex_sc,
                group.by = "subclass",
                label = TRUE) + Seurat::NoLegend()


cortex_sc@meta.data %>% group_by(subclass) %>% summarise(n())


# Get Markers
Seurat::Idents(object = cortex_sc) <- cortex_sc@meta.data$subclass
cluster_markers_all <- Seurat::FindAllMarkers(object = cortex_sc,
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = cluster_markers_all,
        file = "data/spotligh_markers.RDS")


# load spatial
if (! "stxBrain" %in% SeuratData::AvailableData()[, "Dataset"]) {
    # If dataset not downloaded proceed to download it
    SeuratData::InstallData("stxBrain")
}

# 10X Genomics Visium Mouse Brain Dataset
anterior <- SeuratData::LoadData("stxBrain", type = "anterior1")
Seurat::SpatialDimPlot(anterior)



# Run Spotlight ----
set.seed(123)
cluster_markers_all <- readRDS("data/spotligh_test.RDS")
spotlight_ls <- spotlight_deconvolution(
    se_sc = cortex_sc,
    counts_spatial = anterior@assays$Spatial@counts,
    clust_vr = "subclass", # Variable in sc_seu containing the cell-type annotation
    cluster_markers = cluster_markers_all, # Dataframe with the marker genes
    cl_n = 250, # number of cells per cell type to use
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
    mutate(subclass = str_replace_all(subclass, " ", ".")) %>%
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
                 values_to = "predictor") %>%
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

# plot roc lul
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
# 1) cell-cell interaction proportions - correlate on TOP hit proportions

n_rank <- 500
liana_format <- liana_res %>%
    mutate(across(c(source, target), ~str_replace_all(.x, " ", "."))) %>%
    dplyr::select(source, target, ends_with("rank"), -c(mean_rank, median_rank)) %>%
    mutate(aggregate_rank = min_rank(aggregate_rank)) %>%
    pivot_longer(-c(source,target),
                 names_to = "method_name",
                 values_to = "predictor") #rank

# load spotlight deconvolution check
spotlight_ls <- readRDS("data/spotlight_ls.rds")
decon_mtrx <- spotlight_ls[[2]]


corr_deconv_mat <- cor(decon_mtrx) %>%
    reshape2::melt() %>%
    as_tibble() %>%
    dplyr::rename(celltype1 = Var1,
                  celltype2 = Var2,
                  estimate = value)



liana_prop <- liana_format %>%
    filter(predictor < n_rank) %>% # predictor = rank
    group_by(method_name) %>%
    mutate(tot_n = n()) %>%
    group_by(source, target, method_name) %>%
    mutate(n = n()) %>%
    group_by(source, target, method_name) %>%
    mutate(prop = n/tot_n) %>%
    dplyr::select(source, target, method_name, prop) %>%
    distinct() %>%
    group_by(method_name) %>%
    mutate(check_prop = sum(prop))




liana_corr <- liana_prop %>%
    left_join(corr_deconv_mat, by = c("source"="celltype1",
                                    "target"="celltype2"))  %>%
    # FILTER AUTOCRINE
    filter(source!=target) %>%
    group_by(method_name) %>%
    group_nest() %>%
    mutate(correlation = data %>%
               map(function(d) d %>%
                       summarise(correlation = cor(prop,
                                                   estimate,
                                                   method = "pearson")))) %>%
    unnest(correlation)
liana_corr


# Whole vector correlation (not possible)
# Fischer's exact test ----
liana_loc <- liana_format %>%
    left_join(corr_deconv_mat, by = c("source"="celltype1",
                                      "target"="celltype2"))  %>%
    # FILTER AUTOCRINE
    filter(source!=target) %>%
    dplyr::mutate(localisation = case_when(estimate >= 0.2 ~ "colocalized",
                                           estimate < 0.2 ~ "not_colocalized"
                                           )) %>%
    mutate(top_or_not = case_when(predictor <= 1000 ~ "top",
                                  predictor > 1000 ~ "not"
    )) %>% ungroup()



factor1 <- liana_loc %>%
    filter(method_name=="natmi.rank") %>%
    filter(top_or_not == "top") %>% pluck("localisation") %>%
    as.factor() %>%
    recode_factor("colocalized" = TRUE, "not_colocalized" = FALSE)

factor2 <- liana_loc %>%
    filter(method_name=="natmi.rank") %>%
    pluck("localisation") %>%
    as.factor() %>%
    recode_factor("colocalized" = TRUE, "not_colocalized" = FALSE)


fisher.test(factor1, factor2)



xd <- liana_loc %>%
    # filter(method_name=="natmi.rank") %>%
    dplyr::select(method_name, localisation, top_or_not)

entity = "interactions"
resource = "SIGNOR"
var = pathway


ligrec_olap$interactions
enrich2(ligrec_olap$interactions, var1 = pathway, var2 = resource)

liana_loc_filt <- liana_loc %>% filter(predictor <= n_rank)
enrich2(liana_loc_filt, method_name, localisation)



enrich2(xd, "localisation")



ranks_counted <- liana_loc %>%
    # count total (Null)
    group_by(method_name, localisation) %>%
    mutate(total = n())  %>%
    # count in x rank (Alt)
    filter(predictor <= n_rank) %>%
    group_by(method_name, localisation) %>%
    mutate(top = n()) %>%
    dplyr::select(method_name, localisation, total, top) %>%
    distinct()

fet_results <- ranks_counted %>%
    group_by(method_name) %>%
    group_nest(.key = "contigency_tables") %>%
    mutate(fet_res = contigency_tables %>%
               map(function(cont_tab) cont_tab %>% enrich3)) %>%
    unnest(fet_res) %>%
    mutate(
        padj = p.adjust(pval, method = "fdr"),
        enrichment = ifelse(
            odds_ratio < 1,
            -1 / odds_ratio,
            odds_ratio
        ) %>% unname
    ) %>%
    dplyr::select(-contigency_tables)
fet_results

# ^ need to do it over range of ranks on the same graph would show how squidpy
# and cellchat don't change as they simply assign too many interactions /w pval = 0,
# and also could serve as a reference for the other methods




require(tidyverse)
require(rlang)


#' Pairwise enrichment between two factors
#'
#' From a data frame with two categorical variables performs Barnard or
#' Fisher tests between all combinations of the levels of the two variables.
#' In each test the contingency table looks like: in \code{var1} belongs to
#' level x or not vs. in \code{var2} belongs to level y or not. By default
#' Fisher test is used simply because Barnard tests take very very long to
#' compute.
#'
#' @param data A data frame with entities labeled with two categorical
#'     variables \code{var1} and \code{var2}.
#' @param var1 Name of the first categorical variable.
#' @param var2 Name of the second categorical variable. Because in
#'     ligand-receptor data it's always `resource` the default is `resource`.
#' @param p_adj_method Adjustment method for multiple testing p-value
#'     correction (see \code{stats::p.adjust}).
#' @param test_method Characted: either "barnard", "barnard2" or "fisher".
#'     "barnard" uses \code{Barnard::barnard.test} while "barnard2" uses
#'     \code{DescTools::BarnardTest}.
#' @param ... Passed to the function executing the Barnard test.
#'
#' @return A data frame with all possible combinations of the two categorical
#'     variables, p-values, adjusted p-values and odds ratios.
#'
#' @importFrom rlang ensym !! !!! ensym quo_text exec
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr filter group_by group_modify mutate last
#' @importFrom purrr cross_df map
#' @importFrom Barnard barnard.test
#' @importFrom DescTools BarnardTest
#' @importFrom RCurl merge.list
enrich2 <- function(
    data,
    var1,
    var2 = resource,
    p_adj_method = 'fdr',
    test_method = 'fisher',
    ...
){

    var1 <- ensym(var1)
    var2 <- ensym(var2)
    var1_str <- quo_text(var1)
    var2_str <- quo_text(var2)

    t_f <- c('TRUE', 'FALSE')

    data %<>% dplyr::select(!!var1, !!var2)

    data %>%
        {map(names(.), function(x){unique(.[[x]])})} %>%
        setNames(names(data)) %>%
        cross_df() %>%
        filter(!is.na(!!var1) & !is.na(!!var2)) %>%
        group_by(!!var1, !!var2) %>%
        group_modify(
            function(.x, .y){
                f1 <- factor(data[[var1_str]] == .y[[var1_str]][1], levels = t_f)
                f2 <- factor(data[[var2_str]] == .y[[var2_str]][1], levels = t_f)
                result <- fisher.test(f1, f2)
                odds_ratio <- result$estimate
                if(test_method == 'barnard'){
                    sink('NUL')
                    result <- exec(barnard.test, !!!table(f1, f2), ...)
                    sink()
                }else if(test_method == 'barnard2'){
                    param <-
                        list(...) %>%
                        merge.list(list(fixed = NA, method = 'boschloo'))
                    result <- exec(BarnardTest, f1, f2, !!!param)
                }
                tibble(pval = last(result$p.value), odds_ratio = odds_ratio)
            }
        ) %>%
        ungroup() %>%
        mutate(
            padj = p.adjust(pval, method = p_adj_method),
            enrichment = ifelse(
                odds_ratio < 1,
                -1 / odds_ratio,
                odds_ratio
            ) %>% unname
        )

}



enrich3 <- function(df){
    cont_table <- df %>%
        as.data.frame() %>%
        column_to_rownames("localisation") %>%
        t()

    result <- fisher.test(cont_table)
    tibble(pval = last(result$p.value), odds_ratio = result$estimate)
}



# %>%
#     as.data.frame() %>%
#     column_to_rownames("name")



fisher.test(f1, f2)


dat

dat <- data.frame(
    "smoke_no" = c(7, 0),
    "smoke_yes" = c(2, 5),
    row.names = c("Athlete", "Non-athlete"),
    stringsAsFactors = FALSE
)
colnames(dat) <- c("Non-smoker", "Smoker")

dat


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

