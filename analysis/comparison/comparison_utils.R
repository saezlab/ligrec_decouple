#' Function to generate a list in `MethodSpecifics` object forlumation for each
#' method, using the score_specs from liana:::.score_* and a custom score order
#' to deal with permutation method ties
#'
#' @param liana_all_path path to the object generated for a given dataset
#' as obtained via `liana_wrap`
#' @param .score_spec the score specs to be used
#'
#' @returns a list for each method as a MethodSpecifics object
#'
#' @details options for .score_spec are: `liana:::.score_specs`,
#'  `liana:::.score_housekeep`, `.score_comp`
get_spec_list <- function(liana_all_path,
                          .score_spec){

    readRDS(liana_all_path) %>%
        map2(names(.), function(liana_meth, method_name){

            if(is.null(.score_spec()[[method_name]])){
                warning(str_glue(" NO SCORE FOUND FOR {method_name}"))
                return()
            } else{
                message(str_glue("Now reading {method_name}"))
            }

            methods::new("MethodSpecifics",
                         method_name = method_name,
                         method_results = liana_meth,
                         method_scores=list(
                             .score_spec()[[method_name]]@descending_order
                         ) %>% setNames(.score_spec()[[method_name]]@method_score))
        }) %>% compact()
}

#' Helper function to obtain scores used in the comparison
#'
#' @details same as liana:::.score_specs(), but returns uses the alternative
#' metrics to cellchat and squidpy due to ties
.score_comp <- function(){
    sp <- liana:::.score_specs()
    hs <- liana:::.score_housekeep()

    # Squidpy
    sp$squidpy@method_score <- hs$squidpy@method_score
    sp$squidpy@descending_order <- hs$squidpy@descending_order

    # CellChat
    sp$cellchat@method_score <- hs$cellchat@method_score
    sp$cellchat@descending_order <- hs$cellchat@descending_order

    return(sp)
}


#' Get top hits list
#'
#' @param spec_list list of spec objects with ligrec results
#' @param n_ints n of top integers or fractions
#'
#' @inheritDotParams top_enh
#'
#' @return A list of top hits per tool/tool_parameter
#' @export
get_top_hits <- function(spec_list,
                         n_ints,
                         ...){
    map(n_ints, function(.tn){
        names(spec_list) %>%
            map(function(method_name){

                parnams <- names(spec_list[[method_name]]@method_scores)

                map(parnams, function(parm){
                    spec_list[[method_name]]@method_results %>%
                        map(function(res){
                            parm_order <-
                                spec_list[[method_name]]@method_scores[[parm]]

                            res %>%
                                top_enh(n=if_else(parm_order,
                                                .tn,
                                                -.tn),
                                      wt=parm,
                                      ...) %>%
                                as_tibble() %>%
                                distinct() %>%
                                # decomplexify cellchat output
                                {if(method_name=="cellchat") liana:::decomplexify(., columns = c("ligand", "receptor")) %>%
                                        rename(ligand.complex = ligand_complex,
                                               receptor.complex = receptor_complex) %>%
                                        select(ligand, receptor, source, target, everything()) %>%
                                        ungroup()
                                    else .}

                        })
                }) %>%
                    {if(length(spec_list[[method_name]]@method_scores) > 1)
                        setNames(., str_glue("{method_name}_{parnams}"))
                        else setNames(., method_name)
                    }
            }) %>% setNames(names(spec_list)) %>%
            flatten() #%>%
            # setNames(str_replace_all(names(.),
            #                          pattern = "_",
            #                          replacement = "\\."))
    }) %>% setNames(str_glue("top_{n_ints}"))
}


#' S4 Class used to format benchmark output.
#' @name MethodSpecifics-class
#'
#' @field method_name name of the method (e.g. CellChat)
#' @field method_results Named list of method-resource results
#' @field method_scores Named list of the measures provided by the method
#'  and whether they should be interpreted in descending order as the value
#'
#' @exportClass MethodSpecifics
setClass("MethodSpecifics",
         slots=list(method_name="character",
                    method_results = "list",
                    method_scores="list"))


#' Helper Function to handle specific cases for the different Methods
#'
#' @inheritDotParams dplyr::top_n
#' @inheritDotParams dplyr::top_frac
#' @importFrom dplyr top_n
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#'
#' @return Ordered tibble/df as from top_n
top_enh <- function(...,
                    pval_thresh = 0.05,
                    de_thresh = 0.05,
                    sca_thresh = 0.5,
                    top_fun = "top_n"
                    ){

    elipses <- list(...)
    elipses$wt <- sym(elipses$wt)

    # Filter according to:
    if(elipses$wt == "prob"){ # CellChat Probabilities
        elipses[[1]] %<>% filter(pval <= pval_thresh)
    } else if(elipses$wt == "pval"){ # CellChat pval
        elipses[[1]] %<>% filter(pval <= pval_thresh)
    } else if(elipses$wt == "pvalue"){ # Squidpy pvalue
        elipses[[1]] %<>% filter(pvalue <= pval_thresh)
    } else if (elipses$wt == "LRscore"){ # SCA LRscore
        elipses[[1]] %<>% filter(LRscore >= sca_thresh)
    } else if (elipses$wt == "weight_sc"){ # Connectome DE ligrecs
        elipses[[1]] %<>% filter(p_val_adj.rec <= de_thresh) %>%
            filter(p_val_adj.lig <= de_thresh)
    }
    return(do.call(top_fun, elipses))
}




#' Get binary activity frequencies per cell-cell pair
#'
#' @param sig_list List of significant hits per Method-resource combination.
#'  Named list of methods with each element being a named list of resources.
#' @return A tibble of cell pair frequencies, based on binarized activity
#' @export
get_binary_frequencies <- function(sig_list){
    sig_list %>%
        enframe() %>%
        unnest(value) %>%
        mutate(name = map(names(sig_list), # get combined method and resource names
                          function(m_name){
                              map(names(sig_list[[m_name]]),
                                  function(r_name){
                                      str_glue("{m_name}_{r_name}")
                                  })
                          }) %>% unlist()) %>%
        mutate(value = value %>%
                   map(function(res) res %>%
                           unite(source, target, col = "clust_pair") %>%
                           group_by(clust_pair) %>%
                           summarise(cn = n()) %>%
                           mutate(freq = cn / sum(cn)) %>%
                           select(clust_pair, freq)
                   )) %>%
        unnest(value)
}


#' Convert list with MethodSpecifics objects to a dataframe with ranked
#' cell_pair frequencies.
#'
#' @param spec_list list with appropriately populated MethodSpecifics objects
#' @return a tibble with cell_pair frequencies represented by the ranked
#' normalized scores for each method by cell_pair
#' @details See format_rank_frequencies for details
#' @export
get_rank_frequencies <- function(spec_list){
    names(spec_list) %>%
        map(function(method_name){

            parnams <- names(spec_list[[method_name]]@method_scores)

            map(parnams, function(parm){
                spec_list[[method_name]]@method_results %>%
                    map(function(res){
                        res %>%
                            format_rank_frequencies(score_col=sym(parm),
                                                    .desc_order = spec_list[[method_name]]@method_scores[[parm]])
                    })
            }) %>%
                {if(length(spec_list[[method_name]]@method_scores) > 1)
                    setNames(., str_glue("{method_name}_{parnams}"))
                    else setNames(., method_name)
                }
        }) %>% setNames(names(spec_list)) %>%
        flatten() %>%
        setNames(str_replace_all(names(.),
                                 pattern = "_",
                                 replacement = "\\.")) %>%
        reform_rank_frequencies()
}


#' Helper function to get cell pair activity from rank averages
#'
#' @param result Result for a specific Tool-resource combinations
#' @param score_col Score column name provided by the tool
#' @param .desc_order Whether the most significant hits are in descending order,
#' i.e. the highest is the most sig
#' @details Cell pair ranks are averaged, then converted to z-scores, which
#' are multiplied by -1 as we want the lowest average ranks to have the highest
#' z scores
#' @export
format_rank_frequencies <- function(result, score_col, .desc_order = TRUE){
    result %>%
        as_tibble() %>%
        filter(source!=target) %>%
        select(source, target, ligand, receptor, !!score_col) %>%
        filter(!is.nan(!!rlang::sym(score_col))) %>%
        mutate(edge_rank := if_else(rep(.desc_order, nrow(.)),
                                    min_rank(desc(!!rlang::sym(score_col))),
                                    min_rank(!!rlang::sym(score_col))))  %>%
        unite(source, target, col = "clust_pair") %>%
        group_by(clust_pair) %>%
        summarise(avg_rank = (mean(edge_rank))) %>%
        mutate(freq = -1*scale(avg_rank)[, 1])
}



#' Helper function to convert list with all resources ranked to frequencies df
#' @param frequencies_list list with all resources ranked to frequencies df
#' @return rank frequencies df
reform_rank_frequencies <- function(frequencies_list){

    # Combine all results into tool_resource list
    lnames <- map(names(frequencies_list), function(l_name){
        map(names(frequencies_list[[l_name]]), function(r_name){
            str_glue("{l_name}_{r_name}")
        })
    }) %>% unlist()

    freq_df <- frequencies_list %>%
        purrr::flatten() %>%
        setNames(lnames) %>%
        enframe() %>%
        unnest(value)

    return(freq_df)
}






#' Helper Function to get a binary top hits DF (for all method-resource combos)
#' @param sig_list list of significant hits per method-resource combo
#' @export
get_binary_df <- function(sig_list){
    # get method and resource names combined
    lnames <- map(names(sig_list), function(m_name){
        map(names(sig_list[[m_name]]), function(r_name){
            str_glue("{m_name}⊎{r_name}")
        })
    }) %>%
        unlist()

    # get binarized significant hits list (1 for sig per method, 0 if absent)
    binary_df <- sig_list %>%
        purrr::flatten() %>%
        setNames(lnames) %>%
        prepForUpset() %>%
        as_tibble() %>%
        distinct() %>%
        column_to_rownames("interaction")
}



#' Get (Dis)Similarities per Resource/Method Combinations
#' @param sig_list list of top hits
#' @inheritDotParams proxy::simil
#' @import tibble
#' @import purrr
#' @export
simdist_resmet <- function(sig_list,
                           ...){
    binary_df <- get_binary_df(sig_list)
    excl_res <- c(
        "Reshuffled",
        "Default"
        )

    methods <- names(sig_list)
    resources <- names(sig_list[[1]])
    resources <- resources[!(resources %in% excl_res)] # remove these

    # Get Sim/Diss between Methods
    method_sim <- methods %>% map(function(met){
        binary_df %>%
            select(starts_with(met)) %>%
            select(!ends_with(excl_res))
    }) %>% map(function(met_binary)
        get_simil_dist(
            x = t(met_binary),
            ...)) %>%
        setNames(methods)

    # Get Sim/Diss between Resources
    resource_sim <- resources %>% map(function(resource){
        binary_df %>%
            select(ends_with(resource))
    }) %>% map(function(res_binary)
        get_simil_dist(
            x = t(res_binary),
            ...)) %>%
        setNames(resources)

    # This can be extended with other combinations of (dis)similarity lists

    return(list(
        "meth" = method_sim,
        "reso"= resource_sim))
}



#' Get Similarity/Dissimilarity Stats from lists with binary matrices
#' @param ... Any list with with DFs (typically binarized) for which we wish to
#'  to calculate the mean, median, sd, and length.
#' @return a summary tibble
#' @import tibble
#' @export
list_stats <- function(...){
    args <- list(...)

    # combined (i.e. vectorised simdists)
    df_comb <- args %>%
        enframe(value="simdist") %>%
        mutate(name = str_glue("{name}_comb")) %>%
        mutate(simdist = simdist %>%
                   map(function(sim_mat) as.vector(sim_mat) %>% unlist))

    # averaged simdists per each element in the list
    df_mean <- args %>%
        enframe(value="simdist") %>%
        mutate(name = str_glue("{name}_mean")) %>%
        rowwise() %>%
        mutate(simdist = list(simdist %>%
                                  map(function(sim_mat)
                                      mean(sim_mat))
                              %>% unlist)
        ) %>%
        ungroup()

    # bind and calculate averages
    df_stats <- bind_rows(df_comb, df_mean) %>%
        rowwise() %>%
        mutate(mn = mean(simdist),
               med = median(simdist),
               len = length(simdist),
               sd = sd(simdist),
               .min = min(simdist),
               .max = max(simdist)) %>%
        ungroup()

    return(df_stats)
}


#' Creates a path for a figure output
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang !!! exec
#' @export
figure_path_mr <- function(fname,
                           dir_vec,
                           outdir = NULL,
                           ...){

    args <- list(...)
    outdir <- args$outdir
    args$outdir <- NULL
    fname %<>% {exec(sprintf, ., !!!args)}

     outdir %>%
         `%||%`(options('intercell.fig_meth_res')[[1]]) %>%
         `%||%`(file.path("figures","method_resource")) %T>%
         dir.create(showWarnings = FALSE, recursive = TRUE) %>%
         file.path(fname)

}



#' Get Jaccard indexes between specific combination of resources or methods
#' @param sig_list List of significant hits per Method-resource combination.
#'  Named list of methods with each element being a named list of resources.
#' @param methods character vector of method names
#' @param resources character vector of resource names
#'
#' @return Returns a Jaccard index between all combinations of the supplied
#'   resources and methods
#' @export
get_jacc <- function(sig_list, methods, resources){
    get_binary_df(sig_list) %>%
        dplyr::select(ends_with(resources)) %>%
        select(starts_with(methods)) %>%
        filter(rowMeans(.) != 0) %>%
        t() %>%
        get_simil_dist(sim_dist = "simil", "Jaccard")
}




# Seurat Formatting helpers ----


#' Helper function to convert CRC data to sparse Seurat
#' @param counts_loc expression counts location
#' @param meta_loc location of the metadata file
#' @param save_loc save location
sparsify_to_seurat <- function(counts_loc, meta_loc, save_loc){
    counts <- read_delim(counts_loc,
                         delim = "\t") %>%
        as.data.frame() %>%
        column_to_rownames("Index")

    meta <- read_delim(meta_loc,
                       delim = "\t") %>%
        as.data.frame() %>%
        column_to_rownames("Index")

    CreateSeuratObject(counts = Seurat::as.sparse(crc_korean_counts),
                       project = "10X_CRC") %>%
        Seurat::AddMetaData(meta) %>%
        Seurat::NormalizeData() %>%
        Seurat::FindVariableFeatures() %>%
        saveRDS(., save_loc)
}



#' Helper function to get cell number per cell type
#'
#' @param seurat_path Path to Seurat object of interest
#'
#' @return A tibble with Cell_subtype and Cell Number columns
#' @import Seurat dplyr tibble
#' @export
get_cellnum <- function(seurat_path){
    crc_form <- readRDS(seurat_path)
    crc_meta <- crc_form@meta.data

    # Cell Numbers
    crc_meta %>%
        select(Cell_clusters, Cell_subtype) %>%
        group_by(Cell_subtype) %>%
        summarise(cell_occur = n()) %>%
        arrange(Cell_subtype)
}



#' Helper Function to regularize all methods from 0 to 1
#'
#' @param liana_scores liana scores obtained from `get_score_distributions`
#'
#' @details Keep in mind the way that methods are arranged!!!!
regularize_scores <- function(liana_scores,
                              .score_spec = liana:::.score_specs){

    if(identical(.score_spec, liana:::.score_specs)){
        message("LIANA SCORE SPECS!!!!")
        methods_to_revert <- c("CellChat", "CellPhoneDB", "Aggregated Ranks")
    } else{
        message("NOT LIANA SCORE SPECS!!!!")
        methods_to_revert <- c("Aggregated Ranks")
    }

    liana_scores %>%
        # regularize logFC and Connectome to 0-1 (+every negative is 0)
        mutate(score = if_else(method %in% c("LogFC Mean", "Connectome"),
                               if_else(score < 0,
                                       0,
                                       score/max(score)
                               ),
                               score
        )) %>%
        # revert cellchat and squidpy p-values (1 becomes highest)
        rowwise() %>%
        mutate(score = if_else(method %in% methods_to_revert,
                               1 - score,
                               score)) %>%
        ungroup()
}



#' Function to get Activity per Cell Type heatmap
#'
#' @param ct_tibble tibble with ~'activities' per cell type
#' (relative strength/frequncies)
#' > mr(method⊎resource), cell1_source, cell1_target, cell2_source,...
#' @param cap_value Cap cell fraction (prop cell activity) to a given value
#'
#' @return Cell Type Activity Heatmap
#'
#' @import pheatmap tidyverse
#' @inheritDotParams pheatmap::pheatmap
#' @export
get_ct_heatmap <- function(ct_tibble,
                           cap_value = 1,
                           ...){

    # annotation groups (sequential vectors as in heatmap_binary_list)
    method_groups <- ct_tibble %>%
        separate(mr, into = c("method", "resource"), sep = "⊎") %>%
        pull(method)
    resource_groups <- ct_tibble %>%
        separate(mr, into = c("method", "resource"), sep = "⊎") %>%
        pull(resource)

    # data frame with column annotations.
    # with a column for resources and a column for methods
    annotations_df <- data.frame(Resource = resource_groups,
                                 Method = method_groups)  %>%
        mutate(rn = ct_tibble$mr) %>%
        column_to_rownames("rn")

    annotations_row <- data.frame(cell_cat = colnames(ct_tibble)[-1]) %>%
        separate(cell_cat, sep="_", into = c("Cell", "Category"), remove = FALSE) %>%
        column_to_rownames("cell_cat") %>%
        select(Category)

    # List with colors for each annotation.
    mycolors <- list(Method = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(method_groups))),
                     Resource = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(resource_groups))),
                     Category = c("#E41A1C", "#377EB8"))
    names(mycolors$Resource) <- unique(resource_groups)
    names(mycolors$Method) <- unique(method_groups)
    names(mycolors$Category) <- unique(annotations_row$Category)

    lab_rows <- annotations_row %>%
        rownames_to_column("cellname") %>%
        separate(cellname, into = c("cell", "cat"), sep = "_") %>%
        pull(cell)

    cellfraq_heat <- pheatmap::pheatmap(ct_tibble %>%
                                            column_to_rownames("mr") %>%
                                            t(),
                                        annotation_row = annotations_row,
                                        annotation_col = annotations_df,
                                        annotation_colors = mycolors,
                                        display_numbers = FALSE,
                                        silent = FALSE,
                                        show_colnames = FALSE,
                                        color = colorRampPalette(c("darkslategray2",
                                                                   "violetred2"))(20),
                                        fontsize = 30,
                                        drop_levels = TRUE,
                                        cluster_rows = TRUE,
                                        cluster_cols = TRUE,
                                        border_color = NA,
                                        treeheight_row = 0,
                                        treeheight_col = 100,
                                        ...)
}

#' Helper function to get Frequncies of Interactions per Cell Type
#' @param sig_list list of list with top hist (top_frac/top_n)
#' @param cap_value Cap cell fraction (prop cell activity) to a given value
#'
#' @return returns a tibble with method⊎resource frequncies per `cell pair`
get_ct_frequncies <- function(sig_list,
                              cap_value = 1){
    sig_list %>%
        map(function(db){
            db %>%
                enframe(name = "resource", value = "results") %>%
                mutate(results = results %>% map(function(res) res %>%
                                                     select(source, target))) %>%
                unnest(results) %>%
                mutate_at(.vars = c("source", "target"), ~str_replace(., "_", ".")) %>%
                pivot_longer(cols = c(source, target),
                             names_to = "cat",
                             values_to = "cell") %>%
                group_by(resource) %>%
                mutate(total_count = n()) %>%
                ungroup() %>%
                group_by(resource, cat, cell) %>%
                mutate(cell_count = n()) %>%
                group_by(resource, cat, cell) %>%
                mutate(cell_fraq = cell_count/total_count) %>%
                mutate(cell_fraq = if_else(cell_fraq > cap_value,
                                           cap_value,
                                           cell_fraq)) %>%
                distinct() %>%
                unite(cell, cat, col = "cell_cat") %>%
                pivot_wider(id_cols = resource,
                            names_from = cell_cat,
                            values_from = cell_fraq,
                            values_fill = 0)
        }) %>%
        enframe(name = "method", value = "results_resource") %>%
        mutate(method = recode_methods(method)) %>%
        unnest(results_resource) %>%
        unite(method, resource, col = "mr", sep = "⊎") %>%
        mutate_all(~ replace(., is.na(.), 0))
}


#' @title Function to plot distribution densities
#'
#' @param liana_res_specced liana as a specced list (i.e. output of `get_spec_list`)
#' @param hit_prop proportions/fractions of hits to be used
#' @param resource name of the resource to be used
#' @inheritDotParams passed to `get_top_hits` and `liana_aggregate_enh`
#'
#' @returns a dataframe with method scores
get_score_distributions <- function(liana_res_specced,
                                    hit_prop = 1,
                                    resource = "OmniPath",
                                    ...){

    top_frac_lists <- get_top_hits(liana_res_specced,
                                   n_ints=c(hit_prop),
                                   top_fun = "top_frac",
                                   ...)

    # Transpose to resource-method
    liana_resmeth <- top_frac_lists[[str_glue("top_{hit_prop}")]] %>%
        transpose()

    # Format scores
    liana_scores <- liana_resmeth[[resource]] %>%
        map2(names(.), function(met_res, met_name){
            message(met_name)

            met_res %>%
                rename(score = liana:::.score_specs()[[met_name]]@method_score) %>%
                select(source, target, ligand, receptor, score)
        }) %>%
        enframe(name = "method", value = "results") %>%
        unnest(results)


    # liana aggregate rank
    liana_ag_res <- liana_resmeth[[resource]] %>%
        liana_aggregate_enh(...) %>%
        mutate(method = "aggregate_rank") %>%
        select(method, source, target, ligand, receptor, score=aggregate_rank)

    # append aggragate
    liana_scores <- bind_rows(liana_scores, liana_ag_res) %>%
        mutate(method = recode_methods(method))


    return(liana_scores)

}

#' Helper Function to plot Score distributions
#' @param liana_scores liana scores list obtained via `get_score_distributions`
plot_score_distributions <- function(liana_scores){
    # plot
    p <- liana_scores %>%
        ggplot(aes(x=score, color=method, fill=method)) +
        geom_density() +
        theme(text = element_text(size=16)) +
        facet_wrap(~ method, scales='free', nrow = 1) +
        xlab('Method Scores') +
        ylab("Density")
    theme_bw()

    return(p)
}


#' Function to generate relative strenth per cell type (as source and target)
#'
#' @param liana_all_spec liana as a specced list (i.e. output of `get_spec_list`)
#' @inheritDotParams get_score_distributions
#'
#' @returns a tibble with mr (methodUresource), cell1_source, cell2_source, etc
get_ct_strength <- function(liana_all_spec,
                            ...){

    # Get Available Resources
    resources_used <- liana_all_spec[[1]]@method_results %>% names()
    print(resources_used)


    ct_strength <- map(resources_used, function(resource){
        # scores for that resource
        liana_scores <- get_score_distributions(liana_all_spec,
                                                hit_prop = 1,
                                                resource = resource,
                                                ...)
        # Obtain regularized scores
        liana_scores_regularized <- regularize_scores(liana_scores)

        # relative strength per cell pair
        liana_scores_strength <- liana_scores_regularized %>%
            pivot_longer(cols = c(source, target),
                         names_to = "cat",
                         values_to = "cell") %>%
            unite(cell, cat, col = "cell_cat") %>%
            group_by(method) %>%
            mutate(global_score = mean(score)) %>%
            group_by(method, cell_cat) %>%
            mutate(cp_score = mean(score)) %>%
            mutate(cp_strength = cp_score/global_score) %>%
            ungroup() %>%
            select(method, cell_cat, cp_strength) %>%
            distinct() %>%
            pivot_wider(names_from = cell_cat, values_from = cp_strength)

        return(liana_scores_strength)
    }) %>%
        setNames(resources_used) %>%
        enframe(name = "resource") %>%
        unnest(value) %>%
        unite(method, resource, sep ="⊎", col = "mr")
}
