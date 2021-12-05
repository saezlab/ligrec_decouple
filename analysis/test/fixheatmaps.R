ct_tibble <- ct_strength
ph_data <- ct_tibble %>%
    column_to_rownames("mr") %>%
    t()

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
                             Method = method_groups) %>%
    mutate(rn = ct_tibble$mr) %>%
    column_to_rownames("rn")

# data frame with row annotations.
annotations_row <- data.frame(cell_cat = colnames(ct_tibble)[-1]) %>%
    separate(cell_cat, sep="\\^", into = c("Celltype", "Role"), remove = FALSE) %>%
    column_to_rownames("cell_cat") %>%
    select(Role)

# List with colors for each annotation.
mycolors <- list(Method = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(method_groups))),
                 Resource = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(resource_groups))),
                 Role = c("#E41A1C", "#377EB8"))
names(mycolors$Resource) <- unique(resource_groups)
names(mycolors$Method) <- unique(method_groups)
names(mycolors$Role) <- unique(annotations_row$Role)


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

# define legend params
main_title = "Relative\nStrength"
legend_arg_list <- list(title = main_title, #/Frequency
                        title_gp = gpar(fontsize = 24,
                                        fontface = "bold"),
                        grid_height = unit(60, "mm"),
                        grid_width = unit(16, "mm"),
                        legend_height = unit(50, "mm"),
                        size = unit(16, "mm"),
                        labels_gp = gpar(fontsize = 23),
                        pch=32,
                        direction = "vertical")

# define row annotation
left_ann <- rowAnnotation(df = annotations_row,
                          col = list(Role=c("source" = "blue", "target" = "orange")),
                          simple_anno_size = unit(1.4, "cm"),
                          annotation_name_gp = gpar(fontsize = 24),
                          show_legend = TRUE,
                          annotation_legend_param = list(
                              title_gp = gpar(fontsize = 24,
                                              fontface = "bold"),
                              labels_gp = gpar(fontsize = 23),
                              pch=20))



ht <- ComplexHeatmap::Heatmap(
    ph_data,
    col=colorRampPalette(c("darkslategray2",
                           "violetred2"))(10),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    top_annotation = top_ann,
    left_annotation = left_ann,
    heatmap_legend_param = legend_arg_list,
    column_dend_height = unit(5, "cm"),
    row_names_gp = gpar(fontsize = 24),
    row_labels = gsub("\\.", " ", gsub("[\\^].*", "", rownames(ph_data))),
    show_row_names = TRUE,
    show_column_names = FALSE
    )
ht



get_ct_heatmap(ct_strength,
               cap_value = 1,
               main_title = "xx")
