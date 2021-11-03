#' Helper function to prepare list for UpsetPlot
#' @param named_list a named list with sugnificant LR results
#' @return Matrix/DF of 0 and 1 where 1 means that the interaction
#' is present and 0 means that it is not
#' @importFrom tibble tibble
#' @import dplyr
#' @import purrr
#' @export
prepForUpset <- function(named_list){
  map(names(named_list), function(l_name){
      named_list[[l_name]] %>%
      select(1:4) %>%
      unite("interaction", source, target,
            ligand, receptor, sep="_") %>%
      mutate(!!l_name := 1)
  }) %>% reduce(., full_join, by = "interaction") %>%
    mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
    mutate_at(vars(2:ncol(.)), ~ replace(., . != 0, 1)) %>%
    as.data.frame()
}


#' Helper function to save Upset plot
#' @param upset_df prepForUpset output
#' @param dir path to save upset plot
#' @return Null
#' @import UpSetR
#' @export
plotSaveUset <- function(upset_df, file_name, entity_name){
  up <- upset(upset_df, nsets = ncol(upset_df), order.by = "freq",
        point.size = 8, line.size = 3, text.scale	= 3,
        mainbar.y.label = str_glue("Interactions Intersect ({entity_name})"),
        sets.x.label = "Number of Interactions")
  png(file_name, width = 2100, height = 1300)
  print(up)
  dev.off()

  return(up)
}


#' Heatmap looking at binarized top in per method-resource combinations
#'
#' @param sig_list named list of significant hits. Named list of methods with
#'     each element being a named list of resources
#' @inheritDotParams pheatmap::pheatmap
#'
#' @return A pheatmap showing binary overlap between methods and resources
#'
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom purrr map flatten
#' @importFrom magrittr %>%
#' @importFrom stringr str_glue
#' @export
get_BinaryHeat <- function(sig_list,
                        ...){

  heatmap_binary_df <- get_binary_df(sig_list)

  # annotation groups (sequential vectors as in heatmap_binary_df)
  method_groups <- colnames(heatmap_binary_df) %>%
    enframe() %>%
    separate(value, into = c("method", "resource"), sep = "⊎") %>%
    pull(method)
  resource_groups <- colnames(heatmap_binary_df) %>%
    enframe() %>%
    separate(value, into = c("method", "resource"), sep = "⊎") %>%
    pull(resource)

  # data frame with column annotations.
  # with a column for resources and a column for methods
  annotations_df <- data.frame(Resource = resource_groups,
                               Method = method_groups)  %>%
    mutate(rn = colnames(heatmap_binary_df)) %>%
    column_to_rownames("rn")

  # List with colors for each annotation.
  mycolors <- list(Method = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(method_groups))),
                   Resource = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(resource_groups))))
  names(mycolors$Resource) <- unique(resource_groups)
  names(mycolors$Method) <- unique(method_groups)

  binary_heatmap <- pheatmap::pheatmap(heatmap_binary_df,
                                       annotation_col = annotations_df,
                                       annotation_colors = mycolors,
                                       display_numbers = FALSE,
                                       silent = FALSE,
                                       show_rownames = FALSE,
                                       show_colnames = FALSE,
                                       legend_breaks = 0:1,
                                       fontsize = 34,
                                       drop_levels = TRUE,
                                       cluster_rows = FALSE,
                                       cluster_cols = TRUE,
                                       color = c("gray15", "darkslategray2"),
                                       border_color = NA,
                                       clustering_distance_rows = "binary",
                                       clustering_distance_cols = "binary",
                                       treeheight_row = 0,
                                       treeheight_col = 100,
                                       ...
                                       )

  return(binary_heatmap)
}


#' Helper function to Swap Nested Lists
#' @param sig_list named list of significant hits. Named list of methods with
#' each element being a named list of resources
#' @return A list of resources with each element being a named list of methods
#'
#' @details Swap from nested resource lists to nested method lists,
#'  previously the resources were nested by method, this returns the opposite
#' @export
get_swapped_list <- function(sig_list){

  sig_df <- sig_list %>%
    enframe() %>%
    unnest(value) %>%
    mutate(name = map(names(sig_list), # get combined method and resource names
                      function(m_name){
                        map(names(sig_list[[m_name]]),
                            function(r_name){
                              str_glue("{m_name}⊎{r_name}")
                            })
                      }) %>% unlist()) %>%
    separate(name, into = c("method", "resource"), sep = "⊎") %>%
    mutate(value = value %>% setNames(method)) %>%
    group_by(resource)

  # Keep names of resources
  sig_resource_names <- group_keys(sig_df) %>%
    pull(resource)

  sig_list_resource <- sig_df %>%
    group_split() %>%
    map(function(r_list) # get only list values from resource lists
      r_list %>%
        pull(value)) %>%
    setNames(sig_resource_names)
}




#' PCA plot for Cell type/cluster pair frequency for 'significant/top' hits
#' @param freq_df named list of significant hits. Named list of methods with
#' each element being a named list of resources
#' @return A ggplot2 object
#' @import ggfortify ggplot2 RColorBrewer
#' @export
plot_freq_pca <- function(freq_df){
  # format to df with frequencies and
  # Resource and Method columns as factors
  cell_pair_frequency <- freq_df %>%
    pivot_wider(names_from = clust_pair,
                values_from = freq,
                id_cols = name,
                values_fill = 0) %>%
    as.data.frame() %>%
    separate(name, into = c("Method", "Resource"), remove = FALSE, sep="⊎") %>%
    mutate(Method = factor(Method, # prevent ggplot2 from rearranging
                           )) %>%
    mutate(Resource = factor(Resource)) %>%
    column_to_rownames("name")


  # get PCs
  pca_res <- prcomp(cell_pair_frequency[3:ncol(cell_pair_frequency)])

  # frequency plot
  pca_freq <- autoplot(pca_res, data = cell_pair_frequency,
                       colour = "Method", shape = "Resource",
                       size = 6, position = "jitter") +
    scale_color_manual(values=colorRampPalette(brewer.pal(8, "Dark2"))(nlevels(cell_pair_frequency$Method))) +
    theme_bw(base_size = 26) +
    scale_shape_manual(values=1:nlevels(cell_pair_frequency$Resource))

  return(pca_freq)
}





#' Function to get Activity per Cell Type heatmap
#'
#' @param sig_list A resource-tool list with top hits for each combination.
#' @param cap_value Cap cell fraction (prop cell activity) to a given value
#' @return Cell Type Activity Heatmap
#'
#' @import pheatmap tidyverse
#' @inheritDotParams pheatmap::pheatmap
#' @export
get_activecell <- function(sig_list,
                           cap_value = 1,
                           ...){
  # Split by source/target cell
  top_frac <- sig_list %>%
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
    unnest(results_resource) %>%
    unite(method, resource, col = "mr", sep = "⊎") %>%
    mutate_all(~ replace(., is.na(.), 0))

  # annotation groups (sequential vectors as in heatmap_binary_list)
  method_groups <- top_frac %>%
    separate(mr, into = c("method", "resource"), sep = "⊎") %>%
    pull(method)
  resource_groups <- top_frac %>%
    separate(mr, into = c("method", "resource"), sep = "⊎") %>%
    pull(resource)

  # data frame with column annotations.
  # with a column for resources and a column for methods
  annotations_df <- data.frame(Resource = resource_groups,
                               Method = method_groups)  %>%
    mutate(rn = top_frac$mr) %>%
    column_to_rownames("rn")

  annotations_row <- data.frame(cell_cat = colnames(top_frac)[-1]) %>%
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

  cellfraq_heat <- pheatmap::pheatmap(top_frac %>%
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




#' Jaccard Similarities Heatmap Function
#' @param sig_list list of top ranked hits for each method
#' (i.e. top hits obtained via `get_top_hits`)
#' @inheritDotParams get_simil_dist
#' @export
get_simdist_heatmap <- function(sig_list,
                                ...){

  heatmap_binary_df <- get_binary_df(sig_list)
  args <- list(sim_dist = "simil",
               ...)

  simdif_df <- heatmap_binary_df %>%
    t() %>%
    get_simil_dist(.,
                   ...) %>%
    as.matrix()

  # Assign 1 to diagonal
  # Any other distance and similarity works just fine
  # but Jaccard results in NAs in the diagonal,
  # and as.matrix replaces them with 0s..., while class(dist) is immutable
  if(args$method == "Jaccard" && args$sim_dist == "simil"){
    diag(simdif_df) <- 1
  }


  method_groups <- colnames(heatmap_binary_df) %>%
    enframe() %>%
    separate(value, into = c("method", "resource"), sep = "⊎") %>%
    mutate(method=recode_methods(method)) %>%
    pull(method)
  resource_groups <- colnames(heatmap_binary_df) %>%
    enframe() %>%
    separate(value, into = c("method", "resource"), sep = "⊎") %>%
    pull(resource)

  # data frame with column annotations.
  # with a column for resources and a column for methods
  annotations_df <- data.frame(Resource = resource_groups,
                               Method = method_groups)  %>%
    mutate(rn = colnames(heatmap_binary_df)) %>%
    column_to_rownames("rn")

  # List with colors for each annotation.
  mycolors <- list(Method = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(method_groups))),
                   Resource = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(resource_groups))))
  names(mycolors$Resource) <- unique(resource_groups)
  names(mycolors$Method) <- unique(method_groups)

  # define legend params
  legend_arg_list <- list(title = "Jaccard Index",
                          title_gp = gpar(fontsize = 24, fontface = "bold"),
                          grid_height = unit(60, "mm"),
                          grid_width = unit(16, "mm"),
                          legend_height = unit(50, "mm"),
                          size = unit(16, "mm"),
                          labels_gp = gpar(fontsize = 23),
                          pch=32)

  # define top_annotation
  top_ann <- HeatmapAnnotation(df = annotations_df,
                               col = list(Resource=mycolors$Resource,
                                          Method=mycolors$Method),
                               simple_anno_size = unit(1.4, "cm"),
                               annotation_name_gp = gpar(fontsize = 24),
                               show_legend = TRUE,
                               annotation_legend_param = list(title_gp = gpar(fontsize = 24, fontface = "bold"),
                                                              labels_gp = gpar(fontsize = 23),
                                                              pch=20))

  # define row annotation
  left_ann <- rowAnnotation(df = annotations_df,
                            col = list(Resource=mycolors$Resource,
                                       Method=mycolors$Method),
                            simple_anno_size = unit(1.4, "cm"),
                            annotation_name_gp = gpar(fontsize = 24),
                            show_legend = FALSE)

  ht <- ComplexHeatmap::Heatmap(simdif_df,
                                col=colorRampPalette(c("gray15",
                                                       "darkslategray2"))(20),
                                cluster_rows = FALSE,
                                cluster_columns = FALSE,
                                top_annotation = top_ann,
                                left_annotation = left_ann,
                                heatmap_legend_param = legend_arg_list,
                                show_row_names = FALSE,
                                show_column_names = FALSE)

  return(ht)
}



#' Helper Function to get Similarties and Distances from binary dfs
#' @importFrom proxy dist simil
#' @return a similarity/dissimilarity matrix
#' @export
get_simil_dist <- function(sim_dist = "simil", ...){
  do.call(sim_dist,
          list(...))
}


#' @title Recode method names
#' @param dataset - vector /w method names
recode_methods <- function(methods){
  dplyr::recode(methods,
                "squidpy" = "CellPhoneDB",
                "cellchat" = "CellChat",
                "aggregate_rank" = "Aggregated Ranks",
                "logfc" = "LogFC Mean",
                "cytotalk" = "Cytotalk",

                "natmi" = "NATMI",
                "call.natmi" = "NATMI",
                "call_natmi" = "NATMI",

                "call.connectome" = "Connectome",
                "call_connectome" = "Connectome",
                "connectome" = "Connectome",

                "sca" = "SingleCellSignalR",
                "call.sca" = "SingleCellSignalR",
                "call_sca" = "SingleCellSignalR",

                # RNA-ADT correlation baseline
                "RNA-ADT" = "RNA-ADT Baseline"
  )
}



#' @title Function to plot distribution densities
#'
#' @param liana_res_specced liana as a specced list (i.e. output of `get_spec_list`)
#' @param hit_prop proportions/fractions of hits to be used
#' @inheritDotParams passed to `get_top_hits` and `liana_aggregate_enh`
#'
#' @returns ggplot object
plot_score_distributions <- function(liana_res_specced,
                                     hit_prop = 1,
                                     ...){

  top_frac_lists <- get_top_hits(liana_all_spec,
                                 n_ints=c(hit_prop),
                                 top_fun = "top_frac",
                                 ...)

  # Transpose to resource-method
  liana_resmeth <- top_frac_lists[[str_glue("top_{hit_prop}")]] %>%
    transpose()

  # Format scores
  liana_scores <- liana_resmeth$OmniPath %>%
    map2(names(.), function(met_res, met_name){
      message(met_name)

      met_res %>%
        rename(score = liana:::.score_specs()[[met_name]]@method_score) %>%
        select(source, target, ligand, receptor, score)
    }) %>%
    enframe(name = "method", value = "results") %>%
    unnest(results)


  # liana aggregate rank
  liana_ag_res <- liana_resmeth$OmniPath %>%
    liana_aggregate_enh(...) %>%
    mutate(method = "aggregate_rank") %>%
    select(method, source, target, ligand, receptor, score=aggregate_rank)

  # append aggragate
  liana_scores <- bind_rows(liana_scores, liana_ag_res) %>%
    mutate(method = recode_methods(method))

  # plot
  p <- liana_scores %>%
    ggplot(aes(x=score, color=method, fill=method)) +
    geom_density() +
    theme(text = element_text(size=16)) +
    facet_wrap(~ method, scales='free', nrow = 1) +
    xlab('Method Scores') +
    theme_bw()

  return(p)

}


