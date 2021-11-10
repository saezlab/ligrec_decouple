require(liana)
require(tidyverse)
require(magrittr)
require(RColorBrewer)
require(pheatmap)
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
get_ct_heatmap(ct_strength, cap_value = 1)



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
jaccard_per_mr <- simdist_resmet(top_lists[[top_hits_key]],
                                 sim_dist = "simil",
                                 method = "Jaccard")


# Methods jaccard index across different resources
methods_ji <- jaccard_per_mr$meth %>%
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
    unite(method_resource1, method_resource2, col = "combination")


# Boxplot
ggplot(methods_ji,
       aes(x = method,
           y = jacc,
           color = method
           # color = dataset_type,
           # fill = dataset_type
       )) +
    geom_boxplot(alpha = 0.2,
                 outlier.size = 1.5,
                 width = 0.8)  +
    geom_jitter(aes(shape=combination), size = 5, alpha = 0.3, width = 0.15) +
    scale_shape_manual(values = rep(1:20, len = length(unique(methods_ji$combination)))) +
    # facet_grid(~method_name, scales='free_x', space='free', switch="x") +
    theme_bw(base_size = 24) +
    theme(strip.text.x = element_text(angle = 90),
          axis.text.x = element_text(angle = 45, hjust=1)
    ) +
    # guides(color = "none") +
    guides(fill = "none",
           color = "none",
           shape = "none") +
    ylab("Jaccard index") +
    xlab("Methods")



# Resources jaccard index across different methods
resource_ji <- jaccard_per_mr$reso %>%
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
    unite(resource_method1, resource_method2, col = "combination")


ggplot(resource_ji,
       aes(x = resource,
           y = jacc,
           color = resource
           # color = dataset_type,
           # fill = dataset_type
       )) +
    geom_boxplot(alpha = 0.2,
                 outlier.size = 1.5,
                 width = 0.8)  +
    geom_jitter(aes(shape=combination), size = 5, alpha = 0.3, width = 0.15) +
    scale_shape_manual(values = rep(1:20, len = length(unique(resource_ji$combination)))) +
    # facet_grid(~method_name, scales='free_x', space='free', switch="x") +
    theme_bw(base_size = 24) +
    theme(strip.text.x = element_text(angle = 90),
          axis.text.x = element_text(angle = 45, hjust=1)
    ) +
    # guides(color = "none") +
    guides(fill = "none",
           color = "none",
           shape = "none") +
    ylab("Jaccard index") +
    xlab("Resource")
# these become two figures, by method and by resource


# Get Mean of Means of method across resources, and
# Mean of Means of method across resources
jac_summary <- list_stats(meth = jaccard_per_mr$meth,
                          reso = jaccard_per_mr$reso)
jac_summary
# ^ this becomes a table





