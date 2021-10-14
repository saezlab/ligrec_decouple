#' Function to reshape estimate to long format
#' @param estimate_mat diagonal df/matrix with colocalisation estimates
#' @param z_scale whether to z-tranform the estimate column
#'
#' @details converts the diagonal matrix/df to a long format tibble with
#' celltype1, celltype2, and estimate (e.g. correlation/z-scores) columns
reshape_coloc_estimate <- function(estimate_mat, z_scale = FALSE){
    estimate_mat %>%
        reshape2::melt() %>%
        as_tibble() %>%
        dplyr::rename(celltype1 = Var1,
                      celltype2 = Var2,
                      estimate = value) %>%
        # FILTER AUTOCRINE
        filter(celltype1!=celltype2) %>%
        mutate(estimate = ifelse(rep(z_scale, length(estimate)),
                                 scale(estimate),
                                 estimate))
}



#' Function to perform FET on  colocalized vs not_colocalized interactions
#'
#' @param liana_loc liana results formatted with colocalisation in long format
#' with cols:
#' source, target, method_name, predictor (rank), localisation (response)
#' @param arb_thresh a threshold for the co-localization estimate
#'  (i.e. correlation/z-score above or below X)
#'
#' @details What I do is that, I compare (FET) the proportion of colocalized cell types
#' in the top X interactions to the number of colocalized cell types in all interactions.
#' (top vs total)
#' With the assumption that a good method should have a higher proportion of
#' colocalized cell types in its top rank hits.
#'
#' Cont table example:
#'>    localisation       total   top
#'>       <chr>          <int>    <int>
#'>  1 not_colocalized   299868    748
#'>  2 colocalized        20834    159
run_coloc_fet <- function(liana_loc, n_rank){
    # count total vs top colocalized
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

    # contingency tables for colocalized vs not by method
    cont_tabs <- ranks_counted %>%
        group_by(method_name) %>%
        group_nest(.key = "contigency_tables")
    cont_tabs$contigency_tables[[1]]

    # Fischer's exact test on co-localized in top vs total interactions
    fet_results <- cont_tabs %>%
        # FET
        mutate(fet_res = contigency_tables %>%
                   map(function(cont_tab) cont_tab %>% enrich3)) %>%
        # unnest fet results
        dplyr::select(-contigency_tables) %>%
        unnest(fet_res) %>%
        # modify results
        mutate(
            padj = p.adjust(pval, method = "fdr"),
            enrichment = ifelse(
                odds_ratio < 1,
                -1 / odds_ratio,
                odds_ratio
            ) %>% unname
        ) %>%
        arrange(desc(enrichment))
    fet_results

}



#' Function to Tranform LIANA aggregated results
#' @param liana_agg a long tibble with liana_aggregate results
#'
#' @details converts `liana_aggregate` output to long with rank as predictor:
#'    source              target                   method_name     predictor
#     <chr>                <chr>                     <chr>            <dbl>
# Presomitic.mesoderm Anterior.somitic.tissues aggregate_rank          1
liana_agg_to_long <- function(liana_agg){
    liana_agg %>%
        mutate(across(c(source, target), ~str_replace_all(.x, " ", "."))) %>%
        dplyr::select(source, target, ends_with("rank"), -c(mean_rank, median_rank)) %>%
        mutate(aggregate_rank = min_rank(aggregate_rank)) %>%
        pivot_longer(-c(source,target),
                     names_to = "method_name",
                     values_to = "predictor")
}



#' Load10x_Spatial function made to work with matrix rather than h5,
#' adapted from `Seurat::Load10x_Spatial`
#' @param data.dir directory of the data - need to set the wdir to this..
#' @param filedir directory of the matrix files
#' @param project_dir directory of project (used to re-set it...)
#' @inheritParams Seurat::Load10x_Spatial
Load10X_Spatial_enh <- function(data.dir,
                                filedir = "filtered_count_matrix",
                                assay = 'Spatial',
                                slice = 'slice1',
                                filter.matrix = TRUE,
                                to.upper = FALSE,
                                image = NULL,
                                project_dir){
    if (length(x = data.dir) > 1) {
        warning("'Load10X_Spatial' accepts only one 'data.dir'", immediate. = TRUE)
        data.dir <- data.dir[1]
        }
    setwd(file.path(data.dir, filedir)) # SEURAT IS bad dude, doesn't work with the paths...
    data <-  Seurat::ReadMtx(mtx = "matrix.mtx.gz",
                             cells = "barcodes.tsv.gz",
                             features = "features.tsv.gz",
                             feature.column = 1)
    object <- CreateSeuratObject(counts = data, assay = assay)

        if (is.null(x = image)) {
            image <- Read10X_Image(
                image.dir = file.path(data.dir, 'spatial'),
                filter.matrix = filter.matrix
            )
        } else {
            if (!inherits(x = image, what = "VisiumV1"))
                stop("Image must be an object of class 'VisiumV1'.")
        }
        image <- image[Cells(x = object)]
        DefaultAssay(object = image) <- assay
        object[[slice]] <- image

        setwd(project_dir)

        return(object)
    }
