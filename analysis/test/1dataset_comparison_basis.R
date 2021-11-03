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
liana_all <- readRDS("data/output/liana_all_resources.RDS")

# Ranked Scores according to a set of criteria (here by specificy whenever available)
liana_all_spec <- get_spec_list("data/output/liana_all_resources.RDS",
                                #"data/output/crc_res/liana_crc_all.rds", # crc results
                                .score_spec = liana:::.score_specs) # liana:::.score_*

# Top X proportion of hits (according to the ranking specs above)
top_lists <- get_top_hits(liana_all_spec,
                          n_ints=c(0.01, 0.05), # top 1% and top 5%
                          top_fun = "top_frac",
                          pval_thresh = 1,
                          sca_thresh = 0,
                          de_thresh = 0.05)

# Top X of hits
top_lists_n <- get_top_hits(liana_all_spec,
                            n_ints=c(100, 500), # top 100 and 500
                            top_fun = "top_n",
                            pval_thresh = 1,
                            sca_thresh = 0,
                            de_thresh = 0.05)

# I) Score Distributions -----
# This is the one to be used (full methods - no filtering (only DE for connectome)
plot_score_distributions(liana_all_spec,
                         hit_prop = 1,
                         pval_thresh = 1,
                         sca_thresh = 0,
                         de_thresh = 0.05)

# II) Obtain JI Heatmap ----
p <- get_simdist_heatmap(top_lists_n$top_100, # top hits in this case
                         sim_dist = "simil",
                         method = "Jaccard",
                         diag = TRUE,
                         upper = TRUE,
                         cluster_rows = FALSE,
                         cluster_columns = FALSE)
p

# III) Obtain JI Stats ----
# Get Jaccard Stats
jaccard_per_mr <- simdist_resmet(top_lists_n$top_100,
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


# IV) Interaction Frequencies ----
get_activecell(top_lists_n$top_100,
               cap_value = 1
               )


