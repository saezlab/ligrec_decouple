#' Pipeline Function to build eval curves on cytokine activity as ground truth
#'
#' @param seurat_object Seurat object with celltypes
#' @param liana_res liana results (not aggregated)
#' @param cytosig_net cytosig network formatted obyained from `run_cytosig_eval`
#' @param z_scale z-transform NES scores from decoupleR
#' @param expr_prop minimum proportion of gene expression per cell type
#' @param assay Assay to consider
#' @param sum_count_thresh minimum (summed) counts per gene
#' @param NES_thresh NES threshold to be considered as ground truth
#'
#'
#'  @return a tibble with nested roc, prc, and corr for each method
run_cytosig_eval <- function(seurat_object,
                             liana_res,
                             cytosig_net,
                             z_scale,
                             expr_prop,
                             assay,
                             sum_count_thresh,
                             NES_thresh){

    # aggregate liana
    liana_res <- liana_res %>% liana_aggregate()

    # get pseudobulk and filter
    pseudo <- get_pseudobulk(seurat_object,
                             assay = assay,
                             expr_prop = expr_prop,
                             sum_count_thresh = sum_count_thresh)

    message("Pseudobulk Cytokine Enrichment")
    # Cytokine Activity Enrichment
    pseudo_cytosig <- pseudo %>%
        mutate(cytosig_res = logcounts %>%
                   map(function(logc){
                       run_wmean(
                           logc,
                           cytosig_net,
                           .source = "cytokine",
                           .target = "target",
                           .mor = "mor",
                           .likelihood = "weight",
                           times = 1000,
                           seed = 1234,
                           sparse = TRUE,
                           randomize_type = "cols_independently") %>%
                           # keep only norm weighted mean
                           filter(statistic == "norm_wmean") %>%
                           # rename
                           select(cytokine=source,
                                  NES=score,
                                  p_value) %>%
                           # correct p
                           mutate(adj_pvalue = p.adjust(p_value))
                   }))
    saveRDS(pseudo_cytosig, "data/output/cytosig_out/ER_BRCA_cytosig.RDS")

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
        select(celltype, cytosig_res) %>%
        unnest(cytosig_res) %>%
        left_join(alias_tib, by = "cytokine") %>%
        mutate(aliases = if_else(is.na(aliases),
                                 cytokine,
                                 aliases)) %>%
        select(-cytokine, cytokine=aliases) %>%
        select(cytokine, celltype, NES) %>% # exludes p-vals
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
        mutate(response = if_else(NES > NES_thresh,
                                  1,
                                  0)) %>%
        mutate(response = factor(response, levels = c(1, 0))) # first level is the truth

    message("Calculating AUCs")
    # ROC and Correlations
    cytosig_eval %<>%
        group_by(method_name) %>%
        group_nest(.key = "cyto_liana") %>%
        mutate(roc = map(cyto_liana,
                         function(df) calc_curve(df,
                                                 curve="ROC",
                                                 downsampling = FALSE,
                                                 times = 100,
                                                 source_name = "cytokine_in_target",
                                                 auc_only = FALSE)
        )) %>%
        mutate(prc = map(cyto_liana,
                         function(df) calc_curve(df,
                                                 curve="PR",
                                                 downsampling = TRUE,
                                                 times = 1000,
                                                 source_name = "cytokine_in_target",
                                                 auc_only = FALSE))) %>%
        mutate(corr = cyto_liana %>%
                   map(function(df){
                       cor.test(df[["NES"]],
                                df[["predictor"]],
                                method="spearman",
                                exact = FALSE) %>%
                           broom::tidy()
                   }))
}



#' Function to get filtered, processed, and log-transform pseudobulk counts
#'
#'  @param seurat_object Seurat object with celltypes
#'  @param assay Assay to consider
#'  @param expr_prop minimum proportion of gene expression per cell type
#'  @param sum_count_thresh minimum (summed) counts per gene
#'
#'  @return returns a tibble with nested counts and logcounts per celltype
get_pseudobulk <- function(seurat_obect,
                           assay,
                           expr_prop,
                           sum_count_thresh){
    # show celltypes considered
    levels(Idents(seurat_object)) %>%
        map(function(lev) message(lev))

    # Convert Seurat Object to SCE
    sce <- Seurat::as.SingleCellExperiment(seurat_object, assay = assay)
    colLabels(sce) <- SeuratObject::Idents(seurat_object)
    gc()

    # Pseudobulk - sum all counts by celltype (+ gene expression proportions)
    pseudobulk <- scuttle::summarizeAssayByGroup(sce,
                                                 ids = colLabels(sce),
                                                 assay.type = "counts", # raw counts
                                                 statistics = c("sum", "prop"))
    pseudobulk_expr <- pseudobulk@assays@data$sum %>%
        as_tibble(rownames = "gene") %>%
        pivot_longer(-gene, names_to = "celltype", values_to = "sum_count")
    pseudobulk_prop <- pseudobulk@assays@data$prop.detected %>%
        as_tibble(rownames = "gene") %>%
        pivot_longer(-gene, names_to = "celltype", values_to = "prop")

    # Filter according to expression
    pseudo <- pseudobulk_expr %>%
        left_join(pseudobulk_prop, by = c("gene", "celltype")) %>%
        # filter genes not expressed in at least 10% of cells per celltype
        # and only keep those with summed count of at least 10
        filter(prop >= expr_prop) %>%
        filter(sum_count >= sum_count_thresh) %>%
        select(-prop)

    # looks good
    pseudo %>%
        mutate(celltype = as.factor(celltype)) %>%
        ggplot(aes(x=celltype, y = log2(sum_count))) +
        geom_violin(trim=FALSE)

    # Nest by Celltype, format, and normalize
    pseudo %<>%
        group_by(celltype) %>%
        group_nest(.key = "counts") %>%
        # format and log transform
        mutate(logcounts = counts %>%
                   map(function(c) c %>%
                           as.data.frame() %>%
                           column_to_rownames("gene") %>%
                           as.matrix() %>%
                           log2()))
}



#' Function to load and format cytosig signatures
#'
#' @param cytosig_path path to cytosig matrix
#'
#' @details reads cytosig matrix and converts it to long tibble in decoupleR
#' network format
load_cytosig <- function(cytosig_path = "data/input/cytosig/cytosig_signature_centroid.csv"){
    ## read CytoSig
    # models/signatures were obtained from https://github.com/data2intelligence/CytoSig/tree/master/CytoSig
    cyto_signatures <- read.table(cytosig_path,
                                  header = TRUE,
                                  row.names = 1)

    # Format Cytosig
    cytosig_net <- cyto_signatures %>%
        as.data.frame() %>%
        rownames_to_column("target") %>%
        pivot_longer(-target,
                     names_to = "cytokine",
                     values_to = "weight") %>%
        # keep top 500 genes per cytokine
        group_by(cytokine) %>%
        slice_max(n=500, order_by = abs(weight)) %>%
        select(cytokine, target, weight) %>%
        mutate(mor = if_else(weight>0, 1, -1)) %>%
        mutate(weight = abs(weight)) %>%
        mutate(cytokine = if_else(cytokine=="A",
                                  "Activin_A",
                                  cytokine)) %>%
        ungroup()

    return(cytosig_net)
}
