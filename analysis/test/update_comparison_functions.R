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

# Generate New Output
liana_path <- system.file(package = "liana")
testdata <- readRDS(file.path(liana_path , "testdata",
                              "input", "testdata.rds"))

liana_all <- liana_wrap(testdata,
                        resource = liana::show_resources()[-c(1:13)], # [-c(1:2)] all resources except Default
                        method = c('call_natmi', 'call_connectome', 'logfc',
                                   'cellchat', 'call_sca', 'squidpy'),
                        expr_prop=0)
saveRDS(liana_all, "data/output/liana_all_resources.RDS")


# LOAD New Results from different Method-Resource Combinations
liana_all <- readRDS("data/output/liana_all_resources.RDS")


# Format to list of resource-method aggragates
# i.e. obtain a liana_res aggragated for each resource
n_rank = 200

top_aggs <- liana_all %>% purrr::transpose() %>%
    map(function(liana_res){
        liana_res %>%
            liana_aggregate_enh(filt_de_pvals = FALSE, # we don't filter here as too few
                                cap=n_rank) # n_ranks + 1 (+1 to be filtered)
    }) %>%
    map(function(liana_res)
        liana_res %>%
            mutate(across(c(source, target), ~str_replace_all(.x, " ", "."))) %>%
            dplyr::select(source, target, ends_with("rank"), -c(mean_rank, median_rank)) %>%
            mutate(aggregate_rank = min_rank(aggregate_rank)) %>%
            pivot_longer(-c(source,target),
                         names_to = "method_name",
                         values_to = "predictor") %>%
            group_by(method_name))

# Cell-Cell Pair Proportions in top ranks per resource
top_prop <- top_aggs %>%
    map(function(liana_format){
        liana_format %>%
            filter(predictor < n_rank) %>% # predictor = rank
            group_by(method_name) %>%
            mutate(tot_n = n()) %>%
            group_by(source, target, method_name) %>%
            mutate(n = n()) %>%
            group_by(source, target, method_name) %>%
            mutate(prop = n/tot_n) %>%
            dplyr::select(source, target, method_name, prop) %>%
            distinct() %>%
            group_by(method_name) %>%
            mutate(check_prop = sum(prop))
    })

# ^ keep as is for formatting to PCAs, etc


# Try on old results
#### Load OLD Results from different Method-Resource Combinations
# liana_all_old <- list("cellchat" = readRDS("data/output/crc_res/cellchat_results.rds"),
#                       "call_connectome" = readRDS("data/output/crc_res/conn_results.rds"),
#                       "call_italk" = readRDS("data/output/crc_res/italk_results.rds"),
#                       "call_natmi" = readRDS("data/output/crc_res/natmi_results.rds"),
#                       "call_sca" = readRDS("data/output/crc_res/sca_results.rds"),
#                       "squidpy" = readRDS("data/output/crc_res/squidpy_results.rds"))
# saveRDS(liana_all_old, "data/output/crc_res/liana_crc_all.rds")

# to spec
liana_all_spec <- get_spec_list("data/output/crc_res/liana_crc_all.rds",
                                .score_spec = .score_comp) # liana:::.score_*



# Define the numbers of highest interactions that we wish to explore
# and get a list with each threshold as its element
liana_all_spec$call_italk <- NULL # I changed the name of the score

top_lists <- get_top_hits(liana_all_spec,
                          n_ints=c(100,
                                   250,
                                   500,
                                   1000))

#### Combine all binary results into heatmap (top500)
png(filename = figure_path_mr('crc_binheat_top500.png'),
    width = 3000,
    height = 1700)

p <- get_BinaryHeat(top_lists$top_500)
grid::grid.draw(p$gtable)
invisible(dev.off())

#### Activity per Cell type
##### Inferred as the proportion of interaction edges that stem from Source
# Cell clusters or lead to Target Cell clusters in the highest ranked interactions.
png(filename = figure_path_mr('crc_activityheat_top500.png'),
    width = 3000,
    height = 1700)
# Cap to 0.2 (due to high activity in SCA)
p <- get_activecell(top_lists$top_500, cap_value = 0.2)
grid::grid.draw(p$gtable)
invisible(dev.off())


#### Jaccard Index exploration
# Similarity Heatmap (according to Jaccard index)
p <- get_simdist_heatmap(top_lists$top_500,
                         sim_dist = "simil",
                         method = "Jaccard",
                         diag = TRUE,
                         upper = TRUE)

png(filename = figure_path_mr('crc_jaccard_top500.png'),
    width = 3200,
    height = 2800)
grid::grid.draw(p$gtable)
invisible(dev.off())

# Get Jaccard Stats
jaccard_per_mr <- simdist_resmet(top_lists$top_500,
                                 sim_dist = "simil",
                                 method = "Jaccard")
jac <- list_stats(meth = jaccard_per_mr$meth,
                  reso = jaccard_per_mr$reso)

# Here we obtain the mean Jaccard index per resource when using different methods
# and per method when using different resources
jac


# pairwise JI between methods
methods <- c("CellChat", "Connectome", "iTALK", "NATMI", "SCA", "Squidpy")
methods_jacc <- methods %>%
    combn(2) %>%
    t() %>%
    as_tibble() %>%
    unite(c("V1", "V2"), col = "method_combo") %>%
    mutate(methods = str_split(method_combo, "_")) %>%
    mutate(jacc_mat = pmap(list(.x = method_combo,
                                .y = methods),
                           .f = function(.x, .y){
                               list(.x = get_jacc(top_lists$top_500,
                                                  .y,
                                                  as.character(get_lr_resources()))
                               )
                           }
    )
    ) %>%
    unnest(jacc_mat) %>%
    rowwise() %>%
    mutate(jacc_mean = mean(jacc_mat)) %>%
    ungroup() %>%
    arrange(desc(jacc_mean))

# pairwise similarity between each method combination
# i.e. here we don't consider the universe of method-resource combinations
# but rather just the comparison of any two methods when using each resource
methods_jacc


# Jaccard Index using the Same Method with different resources
#  (i.e. similarity among combinations using the same method)
# or in other words how similar are the results from each method with different resources
methods %>% map(function(met){
    get_jacc(top_lists$top_500,
             met,
             as.character(as.character(get_lr_resources()))) %>%
        mean()
}) %>% setNames(methods)



### Supplementary Plots --------------------------------------------------------
#### Binary Jaccard index heatmaps (Supp interactions: 100, 250, and 1000)
png(filename = figure_path_mr('crc_binheat_top100.png'),
    width = 3000,
    height = 1640)
get_BinaryHeat(top_lists$top_100)
invisible(dev.off())

png(filename = figure_path_mr('crc_binheat_top250.png'),
    width = 3000,
    height = 1590)
get_BinaryHeat(top_lists$top_250)
invisible(dev.off())

png(filename = figure_path_mr('crc_binheat_top1000.png'),
    width = 3000,
    height = 1700)
get_BinaryHeat(top_lists$top_1000)
invisible(dev.off())


### Supplementary Note 2 -------------------------------------------------------
#### Conversion Differences

# Jaccard index between the Inbuilt resource for each method
# and the OmniPath Querried counterpart of the same resource
methods_default <- list("CellChatDB",
                        "Ramilowski2015",
                        "connectomeDB2020",
                        "LRdb") %>%
    setNames(gsub("\\_", ".", names(liana_all_spec[-5]))) %>%
    purrr::transpose() %>%
    flatten()

methods_default %>%
    enframe(name = "method", value = "default") %>%
    unnest(default) %>%
    mutate(default = as.character(str_glue("{method}_{default}"))) %>%
    mutate(inbuilt = stringr::str_replace(default, "\\_..*", "_Default")) %>%
    mutate(jacc = pmap_dbl(list(default, inbuilt),
                           .f = function(x, y){
                               get_binary_df(top_lists$top_1000) %>%
                                   select(c(rlang::as_string(x),
                                            rlang::as_string(y))) %>%
                                   t() %>%
                                   get_simil_dist(sim_dist = "simil", "Jaccard")
                           })
    )

# COMPLEX stuff



### Supplementary Note 3 -------------------------------------------------------
#### Cluster Specificity and Method Dissimilarity
# Load method results with non-specific cluster measures
nonspec_list <- get_spec_list("data/output/crc_res/liana_crc_all.rds",
                              .score_spec = liana:::.score_housekeep)

nonspec_top_list <- get_top_hits(nonspec_list,
                                 n_ints=c(500))


#### Cluster-unspecific measures plots
# Interactions Overlap Binary Heatmap
png(filename = figure_path_mr('crc_nonspec_binheat_top500.png'),
    width = 3000,
    height = 1700)

p1 <- get_BinaryHeat(nonspec_top_list$top_500)
grid::grid.draw(p1$gtable)
invisible(dev.off())


# Activity per Cell Type Heatmap
png(filename = figure_path_mr('crc_nonspec_activityheat_top500.png'),
    width = 3000,
    height = 1700)
p2 <- get_activecell(nonspec_top_list$top_500)
grid::grid.draw(p2$gtable)
invisible(dev.off())

# Jaccard Heatmap
p3 <- get_simdist_heatmap(nonspec_top_list$top_500,
                          sim_dist = "simil",
                          method = "Jaccard",
                          diag = TRUE,
                          upper = TRUE)

png(filename = figure_path_mr('crc_nonspec_jaccard_top500.png'),
    width = 3200,
    height = 2800)
grid::grid.draw(p3$gtable)
invisible(dev.off())

# Jaccard Stats
jaccard_nonspec <- simdist_resmet(nonspec_top_list$top_500,
                                  sim_dist = "simil",
                                  method = "Jaccard")
# Mean Jacc Index per Method and Resource
list_stats(meth = jaccard_nonspec$meth,
           reso = jaccard_nonspec$reso)


