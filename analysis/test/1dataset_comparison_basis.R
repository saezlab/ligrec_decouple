require(liana)
require(tidyverse)
require(magrittr)
require(RColorBrewer)
require(ComplexHeatmap)
require(proxy)
require(UpSetR)

source("analysis/comparison/comparison_utils.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

# LOAD New Results from different Method-Resource Combinations
liana_all <- readRDS("data/output/temp/liana_all_resources.RDS")

# Ranked Scores according to a set of criteria (here by specificy whenever available)
liana_all_spec <- get_spec_list("data/output/temp/liana_all_resources.RDS",
                                #"data/output/crc_res/liana_crc_all.rds", # crc results
                                .score_spec = liana:::.score_specs) # liana:::.score_*

# params
top_x <- 0.05
top_fun = "top_frac"
pval_thresh = 1
sca_thresh = 0
de_thresh = 1
resource = "OmniPath"

# in function
top_hits_key <- str_glue({"top_{top_x}"})


# Top X proportion of hits (according to the ranking specs above)
top_lists <- get_top_hits(liana_all_spec,
                          n_ints=top_x, # top 1% and top 5%
                          top_fun = top_fun,
                          pval_thresh = pval_thresh,
                          sca_thresh = sca_thresh,
                          de_thresh = de_thresh)

# Top X of hits
# top_lists_n <- get_top_hits(liana_all_spec,
#                             n_ints=c(100, 500), # top 100 and 500
#                             top_fun = "top_n",
#                             pval_thresh = 1,
#                             sca_thresh = 0,
#                             de_thresh = 0.05)


# I) Score Distributions -----
# obtain Per method list (NOT per resource - too many and they are the same...)
liana_scores <- get_score_distributions(liana_all_spec,
                                        hit_prop = 1,
                                        resource = resource,
                                        pval_thresh = pval_thresh,
                                        sca_thresh = sca_thresh,
                                        de_thresh = de_thresh)
liana_scores
saveRDS(liana_scores, str_glue("data/output/temp/liana_scores_{resource}.RDS"))


# This is the one to be used (full methods - no filtering (only DE for connectome)
# A) Density Distributions
plot_score_distributions(liana_scores)

# II) Interaction Relative Strength per Cell Type -----

# Test Interactions strength per Cell Type and per Cell Pair
# Obtain regularized scores
liana_scores_regularized <- regularize_scores(liana_scores,
                                              .score_spec = liana:::.score_specs)

# relative strength per cell pair
liana_scores_strength <- liana_scores_regularized %>%
    group_by(method) %>%
    mutate(global_score = mean(score)) %>%
    group_by(method, source, target) %>%
    mutate(cp_score = mean(score)) %>%
    mutate(cp_strength = cp_score/global_score) %>%
    ungroup() %>%
    select(method, source, target, cp_strength) %>%
    unite(source, target, col = "cell_pair") %>%
    distinct() %>%
    arrange(desc(cp_strength)) %>%
    pivot_wider(names_from = cell_pair,
                values_from = cp_strength) %>%
    as.data.frame() %>%
    column_to_rownames("method") %>%
    as.matrix()

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
    pivot_wider(names_from = cell_cat, values_from = cp_strength) %>%
    as.data.frame() %>%
    column_to_rownames("method") %>%
    as.matrix()


# ^ we get some NAs in the test data (Connectome no sig after filtering)

# Relative regularized strength per cell type pair
# The idea is to show that different methods assign different
# interaction strength to different cell type pairs
# to be split to source and target (as in the other cell pair frequency heatmap)
# i.e. I should merge them
ComplexHeatmap::Heatmap(t(liana_scores_strength))
# get cp_strength
pval_thresh = 1
sca_thresh = 0
de_thresh = 0.05

ct_strength <- get_ct_strength(liana_all_spec,
                               pval_thresh = pval_thresh,
                               sca_thresh = sca_thresh,
                               de_thresh = de_thresh)
get_ct_heatmap(ct_strength, cap_value = 1, main_title="Relative\nStrength")


ct_tibble <- ct_strength
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

annotations_row <- data.frame(cell_cat = colnames(ct_tibble)[-1]) %>%
    separate(cell_cat, sep="\\^", into = c("Cell", "Category"), remove = FALSE) %>%
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

ph_data <- ct_tibble %>%
    column_to_rownames("mr") %>%
    t()


# names(ct_tibble)[-1] <- str_to_title(gsub("[\\^].*", "", names(ct_tibble)[-1]))
cellfraq_heat <- pheatmap::pheatmap(ph_data,
                                    annotation_row = annotations_row,
                                    annotation_col = annotations_df,
                                    annotation_colors = mycolors,
                                    display_numbers = FALSE,
                                    silent = FALSE,
                                    show_colnames = FALSE,
                                    show_rownames = TRUE,
                                    color = colorRampPalette(c("darkslategray2",
                                                               "violetred2"))(20),
                                    fontsize = 26,
                                    drop_levels = TRUE,
                                    cluster_rows = FALSE,
                                    cluster_cols = TRUE,
                                    border_color = NA,
                                    labels_row = gsub("[\\^].*", "", rownames(ph_data)[-1]),
                                    treeheight_row = 0,
                                    treeheight_col = 100
)



# III) Interaction Frequencies per Cell Type ----
# This works TOP x
ct_frequncies <- get_ct_frequncies(top_lists[[top_hits_key]])

get_ct_heatmap(ct_frequncies, cap_value = 1)



# IV) Obtain JI Heatmap ----
p <- get_simdist_heatmap(top_lists[[top_hits_key]], # top hits in this case
                         sim_dist = "simil",
                         method = "Jaccard",
                         diag = TRUE,
                         upper = TRUE,
                         cluster_rows = TRUE,
                         cluster_columns = TRUE)
p

# V) Obtain JI Stats ----
# Get Jaccard Stats
binary_df <- get_binary_df(top_lists[[top_hits_key]])
jaccard_per_mr <- simdist_resmet(sig_list = top_lists[[top_hits_key]],
                                 binary_df = binary_df,
                                 sim_dist = "simil",
                                 method = "Jaccard")



across_methods_ji <- jaccard_per_mr$reso %>%
    compact() %>%
    map2(names(.), function(reso_ji, reso_name){
        reso_ji %>%
            as.matrix() %>%
            as.data.frame() %>%
            as_tibble(rownames="resource_method1") %>%
            pivot_longer(-resource_method1,
                         names_to = "resource_method2",
                         values_to = "jacc") %>%
            distinct() %>%
            filter(resource_method1!=resource_method2) %>%
            mutate(resource = reso_name)
    }) %>%
    bind_rows() %>%
    unite(resource_method1, resource_method2, col = "combination") %>%
    mutate(resource = recode_resources(resource))
methods_jaccbox <- jacc_1d_boxplot(across_methods_ji, entity="resource")

across_resources_ji <- jaccard_per_mr$meth %>%
    compact() %>%
    map2(names(.), function(met_ji, met_name){
        met_ji %>%
            as.matrix() %>%
            as.data.frame() %>%
            as_tibble(rownames="method_resource1") %>%
            pivot_longer(-method_resource1,
                         names_to = "method_resource2",
                         values_to = "jacc") %>%
            distinct() %>%
            filter(method_resource1!=method_resource2) %>%
            mutate(method = met_name)
    }) %>%
    bind_rows() %>%
    unite(method_resource1, method_resource2, col = "combination") %>%
    mutate(method = recode_methods(method))

resources_jaccbox <- jacc_1d_boxplot(across_resources_ji, entity="method")



# Get Mean of Means of method across resources, and
# Mean of Means of method across resources
jac_summary <- list_stats(meth = jaccard_per_mr$meth,
                          reso = jaccard_per_mr$reso)
jac_summary
# ^ this becomes a table





