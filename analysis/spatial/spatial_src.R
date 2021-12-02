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

    # Fischer's exact test on co-localized in top vs total interactions
    fet_results <- cont_tabs %>%
        # FET
        mutate(fet_res = contigency_tables %>%
                   map(function(cont_tab){
                       # handle case where no colocolalized are in top X ranks
                       if(nrow(cont_tab) == 1){
                           if(cont_tab %>% pluck("localisation")=="not_colocalized"){
                               message("Only NOT colocalized are present")
                               tibble(pval=1,
                                      odds_ratio=-9999)
                           }
                       } else{ # run enrichment
                           cont_tab %>% enrich3
                       }
                   })) %>%
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



#' Helper function to check localisation by dataset
#' @param lr_coloc prepped tibble as from analyse_spatial section
check_coloc <- function(lr_coloc){
    lr_coloc %>%
        group_by(dataset, localisation) %>%
        summarise(count_loc = n()) %>%
        pivot_wider(names_from = localisation,
                    values_from = count_loc) %>%
        mutate(proportion = colocalized/not_colocalized)
}


#' Helper Function to run FET and fortmat to boxplot data
#' @param lr_coloc ligand-receptor colocalized data
#' @param n_ranks vector with ranks c(50,100,250,etc)
#' @returns tibble used for get_boxplot
get_fet_boxplot_data <- function(lr_coloc, n_ranks){

    # RUN FET and format
    fet_res <- lr_coloc %>%
        group_by(dataset) %>%
        group_nest(.key="lr_coloc") %>%
        mutate(fet = lr_coloc %>%
                   map(function(liana_loc){
                       fets <- n_ranks %>%
                           map(function(n_rank){
                               run_coloc_fet(liana_loc = liana_loc,
                                             n_rank = n_rank) %>%
                                   mutate(n_rank = n_rank)
                           }) %>%
                           bind_rows()
                   })) %>%
        select(-lr_coloc) %>%
        unnest(fet)

    # Address ties issue
    ties_issue <- lr_coloc %>%
        group_by(method_name, dataset) %>%
        filter(predictor < max(n_ranks)) %>%
        summarise(count_ranks = n()) %>%
        filter(count_ranks > max(n_ranks))

    # Fix issue with single-class cases
    n_issues <- fet_res %>%
        filter(odds_ratio==-9999)

    message(str_glue("{nrow(n_issues)} of 1s appended due to single-class cases"))
    fets <- fet_res %>%
        mutate(odds_ratio = if_else(odds_ratio==-9999,
                                    1,
                                    odds_ratio))

    # Format for boxplot
    boxplot_data <- fets %>%
        # replace squidpy_cellchat ranks
        left_join(ties_issue) %>%
        mutate(n_rank = ifelse(is.na(count_ranks),
                               n_rank,
                               count_ranks
        )) %>%
        select(-count_ranks) %>%
        mutate(n_rank = as.factor(n_rank)) %>%
        distinct_at(.vars = c("method_name", "enrichment", "n_rank", "dataset"),
                    .keep_all = TRUE) %>%
        # recode methods
        mutate(method_name = gsub("\\..*","", method_name)) %>%
        mutate(method_name = recode_methods(method_name)) %>%
        # recode datasets
        mutate(dataset = recode_datasets(dataset))

    return(boxplot_data)
}


#' Function to get boxplot for FET odds ratios
#' @param boxplot_data tibble in the following format:
#' > method_name            pval odds_ratio padj  enrichment n_rank    dataset
#' >     <chr>             <dbl>   <dbl>    <dbl>   <dbl>     <fct>     <chr>
#' > Aggregated Ranks    2.17e- 9  9.21  5.07e- 9   9.21       50  Cortex Anterior 1
#' @returns a ggplot odds-ratio boxplot across methods and datasets
get_spatial_boxplot <- function(boxplot_data){

    if(length(unique(boxplot_data$dataset)) > 3){
       box_or_not <- geom_boxplot(alpha = 0.15,
                     outlier.size = 1.5,
                     width = 0.2,
                     show.legend = FALSE)
    } else{
        box_or_not <- NULL
    }

    # plot Enrichment of colocalized in top vs total
    boxplot <- ggplot(boxplot_data,
                      aes(x = n_rank, y = odds_ratio,
                          color = method_name)) +
        box_or_not +
        geom_jitter(aes(shape = dataset), width = 0, size = 5) +
        facet_grid(~method_name, scales='free_x', space='free', switch="x") +
        theme_bw(base_size = 24) +
        geom_hline(yintercept = 1, colour = "pink",
                   linetype = 2, size = 1.5) +
        theme(strip.text.x = element_text(angle = 90),
              axis.text.x = element_text(angle = 90, hjust=1)
        ) +
        labs(shape=guide_legend(title="Dataset")) +
        ylab("Odds Ratio") +
        xlab("#Ranks Considered") +
        guides(color = "none")

    return(boxplot)
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


#' Prepare brca atlases and markers
#' @param slide_subtype ER or TNBC
#' @param cluster_key cluster key..
#'
#' @details format, preprocess, and get markers for BRCA type atlases
prep_brca_atlases <- function(slide_subtype, cluster_key = "celltype_major"){

    # Define deconvolution results and input dir
    deconv_directory <- file.path(project_dir, brca_dir,
                                  "deconv", str_glue("{slide_subtype}_{cluster_key}"))
    dir.create(deconv_directory, showWarnings = FALSE)
    message(str_glue("Created: {deconv_directory}"))

    message(file.path(project_dir, brca_dir,
                      str_glue("brca_{slide_subtype}_seurat.RDS")))

    # Load BRCA Subtype Atlas
    atlas_object <- readRDS(file.path(project_dir, brca_dir,
                                      str_glue("brca_{slide_subtype}_seurat.RDS")))

    # subsample
    message("Formatting metadata")
    submeta <- atlas_object@meta.data %>%
        rownames_to_column("barcode") %>%
        group_by(!!sym(cluster_key)) %>%
        mutate(ncels_by_group = n()) %>%
        filter(ncels_by_group >= 25) %>% # as in Wu et al., 2021
        as.data.frame() %>%
        mutate({{ cluster_key }} := as.factor(as.character(.data[[cluster_key]]))) %>%
        as.data.frame() %>%
        column_to_rownames("barcode") %>%
        ungroup()

    # reassign meta/clusters
    atlas_object@meta.data <- submeta
    Idents(atlas_object) <- submeta[[cluster_key]]
    atlas_object[[cluster_key]] <- Idents(atlas_object)
    # subset to meta
    atlas_object <- subset(atlas_object, cells = rownames(submeta))

    # Normalize object
    set.seed(1234)
    message("Normalizing")
    atlas_object %<>%
        Seurat::SCTransform(verbose = FALSE,
                            conserve.memory = TRUE) %>%
        Seurat::RunPCA(verbose = FALSE) %>%
        Seurat::RunUMAP(dims = 1:30, verbose = FALSE)

    # save the object
    seurat_path <- file.path(deconv_directory,
                             str_glue("{slide_subtype}_{cluster_key}_seurat.RDS"))
    saveRDS(atlas_object, seurat_path)

    # Find Markers
    message("Running find markers")
    cluster_markers_all <-
        Seurat::FindAllMarkers(object = atlas_object,
                               assay = "SCT",
                               slot = "data",
                               verbose = TRUE,
                               only.pos = TRUE)
    # save markers
    markers_path <- file.path(deconv_directory,
                              str_glue("{slide_subtype}_{cluster_key}_markers.RDS"))
    saveRDS(cluster_markers_all, markers_path)

    return()
}





#' Function to loop over slide_name and slide_subtype
#' @param slide_name name of the visium slide
#' @param slide_subtype cancer subtype of the visium slide
#' @param cluster_key name of the cluster to be used for deconvoluton
#' @param n_cells number of cells to be used
#' @param subsampled whether to run with subsampled single-cell atlas
deconv_brca_slides <- function(slide_name,
                               slide_subtype,
                               cluster_key = "celltype_minor",
                               n_cells = 100,
                               subsampled = TRUE){

    # deconvolution subdirectory
    deconv_directory <- file.path(project_dir, brca_dir,
                                  "deconv", str_glue("{slide_subtype}_{cluster_key}"))
    message(str_glue("Now running {slide_name}({slide_subtype}) using {cluster_key} with {n_cells} cells"))
    message(str_glue("To be Saved/Loaded in: {deconv_directory}"))

    # Load Subsampled atlas and corresponding markers
    atlas_object <-
        readRDS(file.path(deconv_directory,
                          str_glue("{slide_subtype}_{cluster_key}_seurat.RDS")))
    cluster_markers_all <-
        readRDS(file.path(deconv_directory,
                          str_glue("{slide_subtype}_{cluster_key}_markers.RDS")))
    deconv_res_path <-
        file.path(deconv_directory,
                  str_glue("{slide_name}_{slide_subtype}_{cluster_key}_deconv.RDS"))

    # load visium object
    vis_obj <- Load10X_Spatial_enh(file.path(project_dir,
                                             brca_dir,
                                             "/visium_rearranged",
                                             slide_name),
                                   project_dir = project_dir)

    # Run deconvolution as from tutorial and paper
    set.seed(1234)
    spotlight_ls <- spotlight_deconvolution(
        se_sc = atlas_object,
        counts_spatial = vis_obj@assays$Spatial@counts,
        clust_vr = cluster_key, # cell-type annotations
        cluster_markers = cluster_markers_all, # df with marker genes
        cl_n = n_cells, # number of cells per cell type
        hvg = 3000, # Number of HVG
        ntop = NULL, # How many of the marker genes to use (by default all)
        transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorization and NLS
        method = "nsNMF",
        assay = "SCT"
    )

    saveRDS(spotlight_ls,
            deconv_res_path)

    return()
}


#' Plot Spatial Fet Results
#' @param lr_coloc_path path to FET enrichments of co-localized or not
#'
#' @returns a ggplot2 object
get_spatial_bigbox <- function(lr_coloc_path){
    type_lr_coloc <- readRDS(lr_coloc_path)
    ggplot(type_lr_coloc,
           aes(x = n_rank,
               y = odds_ratio,
               color = dataset_type,
               fill = dataset_type)) +
        geom_boxplot(alpha = 0.2,
                     outlier.size = 1.5,
                     width = 0.8)  +
        # geom_point(aes(shape = dataset), size = 2) +
        scale_shape_manual(values = rep(1:12, len =  length(unique(type_lr_coloc$dataset)))) +
        facet_grid(~method_name, scales='free_x', space='free', switch="x") +
        theme_bw(base_size = 32) +
        geom_hline(yintercept = 1, colour = "black",
                   linetype = 2, size = 1.3) +
        theme(strip.text.x = element_text(angle = 90),
              axis.text.x = element_text(angle = 90,
                                         hjust=1, vjust = 0.5,
                                         size = 20)
        ) +
        # guides(color = "none") +
        labs(shape = guide_legend(title="Visium slide"),
             color = guide_legend(title="Dataset type")) +
        guides(fill = "none") +
        ylab("Odds Ratio") +
        xlab("#Ranks Considered")
}
