#' Helper function to prepare list for UpsetPlot
#' @param named_list a named list with significant LR results
#' @return Matrix/DF of 0 and 1 where 1 means that the interaction
#' is present and 0 means that it is not
#' @importFrom tibble tibble
#' @import dplyr
#' @import purrr
#' @export
prepForUpset <- function(named_list){
  map(names(named_list), function(l_name){

    if(nrow(named_list[[l_name]])==0){
      message(str_glue("{l_name} has 0 hits"))
      return(tibble(interaction=character()))
    }

      named_list[[l_name]] %>%
        select(source, target, ligand, receptor) %>%
        unite("interaction", source, target,
              ligand, receptor, sep="_") %>%
        mutate(!!l_name := 1) %>%
        distinct() %>%
        # Convert to DT
        data.table::setDT() %>%
        data.table::setkey(., interaction)
  }) %>%
    Reduce(function(...) merge(..., all = TRUE), .) %>%
    as.data.frame() %>%
    mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
    mutate_at(vars(2:ncol(.)), ~ replace(., . != 0, 1))
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



#' Jaccard Similarities Heatmap Function
#' @param sig_list list of top ranked hits for each method
#' (i.e. top hits obtained via `get_top_hits`)
#' @inheritDotParams get_simil_dist
#' @export
#' Jaccard Similarities Heatmap Function
#'
#' @param sig_list list of top ranked hits for each method
#' (i.e. top hits obtained via `get_top_hits`)
#'
#' @inheritDotParams get_simil_dist
#' @export
get_simdist_heatmap <- function(sig_list,
                                binary_df = NULL,
                                simdif_df = NULL,
                                cluster_rows = FALSE,
                                cluster_columns = FALSE,
                                ...){


  # This allows binary_df to be passed, if already calculated
  if(is.null(binary_df)){
    heatmap_binary_df <- get_binary_df(sig_list)
  } else{
    heatmap_binary_df <- binary_df
  }

  # This if is used to skip the calculations below
  # I use it for the mean JI Heatmap across datasets
  if(is.null(simdif_df)){
    # Calculate data from binary df
    simdif_df <- get_heatmap_data(heatmap_binary_df,
                                  ...)
  }

  method_groups <- colnames(heatmap_binary_df) %>%
    enframe() %>%
    separate(value, into = c("method", "resource"), sep = "⊎") %>%
    mutate(method = recode_methods(method)) %>%
    pull(method)
  resource_groups <- colnames(heatmap_binary_df) %>%
    enframe() %>%
    separate(value, into = c("method", "resource"), sep = "⊎") %>%
    mutate(resource = recode_resources(resource)) %>%
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
  legend_arg_list <- list(title = "Jaccard\nIndex",
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
                               annotation_legend_param = list(title_gp = gpar(fontsize = 24,
                                                                              fontface = "bold"),
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
                                cluster_rows = cluster_rows,
                                cluster_columns = cluster_columns,
                                top_annotation = top_ann,
                                left_annotation = left_ann,
                                heatmap_legend_param = legend_arg_list,
                                column_dend_height = unit(5, "cm"),
                                show_row_names = FALSE,
                                show_column_names = FALSE)

  return(ht)
}



#' Helper function to estimate and format JI heatmap data
get_heatmap_data <- function(heatmap_binary_df,
                             ...){
  simdif_df <- heatmap_binary_df %>%
    t() %>%
    get_simil_dist(.,
                   ...) %>%
    as.matrix()

  # Assign 1 to diagonal
  # Any other distance and similarity works just fine
  # but Jaccard results in NAs in the diagonal,
  # and as.matrix replaces them with 0s..., while class(dist) is immutable
  args <- list(sim_dist = "simil",
               ...)
  if(args$method == "Jaccard" && args$sim_dist == "simil"){
    diag(simdif_df) <- 1
  }

  return(simdif_df)
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
#' @param methods - vector /w method names
recode_methods <- function(methods){
  dplyr::recode(methods,
                "squidpy" = "CellPhoneDB",
                "cellphonedb" = "CellPhoneDB",
                "cellchat" = "CellChat",
                "aggregate_rank" = "Consensus*",
                "Aggregated Ranks" = "Consensus*",
                "logfc" = "LogFC Mean",
                "cytotalk" = "Crosstalk scores",

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



#' Recode colours for Eval figs
recode_colours <- function(colours){
  dplyr::recode(colours %>% sort(),
                "ER+ Breast Cancer" = "#FFA44FFA",
                "Triple N. Breast Cancer" = "#7264B9",
                "HER2+ Breast Cancer" =  "#DEDA00",
                "Brain Cortex" = "#1B9E77",
                "seqFISH" = "#73477A",
                "merFISH" = "#BD6A40"
  )
}

#' @title Recode dataset names
#' @param dataset - vector /w dataset names
recode_datasets <- function(datasets){
  dplyr::recode(datasets,
                # CITE-Seq
                "10k_malt" = "10kMALT",
                "10k_pbmcs" = "10kPBMCs",
                "5k_pbmcs" = "5kPBMCs ",
                "5k_pbmcs_nextgem" = "5kPBMCs (nextgem)",
                "cmbcs" = "3kCBMCs",
                "spleen_lymph_101" = "SLN111",
                "spleen_lymph_206" = "SLN208",

                # Mouse Brain Cortex visium
                "anterior1" = "Cortex Anterior 1",
                "anterior2" = "Cortex Anterior 2",
                "posterior1" = "Cortex Posterior 1",
                "posterior2" = "Cortex Posterior 2",

                "1142243F" = "TNBC1 (1142243F)",
                "1160920F" = "TNBC2 (1160920F)",
                "CID4290" = "ER1 (CID4290)",
                "CID4465" = "TNBC3 (CID4465)",
                "CID4535" = "ER2 (CID4535)",
                "CID44971" = "TNBC4 (CID44971)",
                "ER" = "ER+ Breast Cancer",
                "ER+ BRCA" = "ER+ Breast Cancer",
                "TNBC BRCA" = "Triple N. Breast Cancer",
                "TNBC" = "Triple N. Breast Cancer",
                "HER2" = "HER2+ Breast Cancer",
                "HER2+ BRCA" = "HER2+ Breast Cancer"
  )
}



# Function to Combine the plots
cyto_space_patch <- function(cytosig_p,
                             space_p,
                             path
){
  cairo_pdf(path,
            height = 24,
            width = 22,
            family = 'DINPro')
  print((cytosig_p / space_p) +
          plot_layout(guides = 'keep', heights = c(1, 1)) +
          plot_annotation(tag_levels = 'A',
                          tag_suffix = ')') &
          theme(plot.tag = element_text(face = 'bold',
                                        size = 40)))
  dev.off()
}


