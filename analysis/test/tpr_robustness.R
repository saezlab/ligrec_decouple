require(tidyverse)
source("analysis/robustness/Code/Utilities/Iterator_Functions.R")
source("analysis/robustness/Code/Utilities/User_Outputs_and_Plots.R")
source("src/plot_utils.R")


## Convert RData to RDS
tpr_tobe_calculated <-
    list(
        # Subset cells
        "Cells Removed (%)" = "analysis/robustness/Outputs/Cluster_Reshuffling/Iterator_Results_CR_seurat_pbmc_subset_top250_2021-12-10_14-55.RData",
        # Reshuffle Cells
        "Reshuffled Cell Labels (%)" = "analysis/robustness/Outputs/Cluster_Reshuffling/Iterator_Results_CR_seurat_pbmc_reshuffle_top250_2021-12-10_14-55.RData",
        # Shuffle Resource (preserve top 500)
        "Interactions Replaced (%)" =  "analysis/robustness/Outputs/Resource_Dilution/Iterator_Results_RD_seurat_pbmc_rand_topo_variable_top250_2021-12-10_14-55.RData",
        # Shuffle Resource (don't preserve/indiscriminant)
         "Indiscriminantly Replaced Interactions (%)" = "analysis/robustness/Outputs/Resource_Dilution/Iterator_Results_RD_seurat_pbmc_mod_baseline_rand_topo_variable_top250_2021-12-10_14-55.RData"
        )


# Subset
result_path <- "analysis/robustness/Outputs/Cluster_Reshuffling/Iterator_Results_CR_seurat_pbmc_subset_top250_2021-12-10_14-55.RData"
load(result_path)
calculate_tpr(result_path, analysis_focus = "cluster") %>%
    plot_tpr(descript = "Cells Removed (%)")


# Reshuffle
result_path <- "analysis/robustness/Outputs/Cluster_Reshuffling/Iterator_Results_CR_seurat_pbmc_reshuffle_top250_2021-12-10_14-55.RData"
load(result_path)
calculate_tpr(result_path, analysis_focus = "cluster") %>%
    plot_tpr(descript = "Reshuffled Cell Labels (%)")


# discriminant resource dilution
result_path <- tpr_tobe_calculated$`Interactions Replaced (%)`
load(result_path)
calculate_tpr(result_path, analysis_focus = "resource") %>%
    plot_tpr(descript = "Interactions Replaced (%)")



"analysis/robustness/Outputs/Resource_Dilution/"

# Indiscriminant resource dilution
result_path <- "analysis/robustness/Outputs/Resource_Dilution/Iterator_Results_RD_seurat_pbmc_mod_baseline_rand_topo_variable_top250_2021-12-10_14-55.RData"
load(result_path)
calculate_tpr(result_path, analysis_focus = "resource") %>%
    plot_tpr(descript = "Indiscriminantly Replaced Interactions (%)")


# Generate TPR
load(result_path)




### Combine all results
res_to_tpr <- imap(top_ranks, function(method_res, method_name){

    gt_res <- top_ranks[[method_name]]$Reshuffle_0$Seed_1

    method_res %>%
        enframe(name='shuffle', value = 'results') %>%
        unnest(results) %>%
        mutate(jacc = map_dbl(results, function(res){
            rank_overlap(comparison_ranks = res, main_ranks = gt_res)
        })
        ) %>%
        mutate(tpr = map_dbl(results, function(res){
            tpr_overlap(comparison_ranks = res, main_ranks = gt_res)
        })
        ) %>%
        select(-results)
    }) %>%
    enframe(name = 'Method') %>%
    unnest(value) %>%
    separate(shuffle, into = c("x", "value"), remove = FALSE) %>%
    select(-x) %>%
    mutate(value = as.numeric(value)) %>%
    mutate(Method = recode_methods(Method))


### TPR instead of Jaccard





# TWO LINES
tpr_plot <- convert_to_tpr %>%
    ggplot(aes(x = value,
               y = tpr,
               group = Method))  +
    stat_summary(fun = mean, geom = "line",
                 size = 1.5, color = 'darkblue') +
    geom_errorbar(aes(ymin=tpr-tpr_sd, ymax=tpr+tpr_sd), width=.2,
                  position=position_dodge(0.05)) +
    stat_summary(mapping = aes(x=value , y = fpr),
                 fun = mean, geom = "line",
                 size = 1.5, color = 'darkred') +
    # make sure the scale is the same on both axes, since both are in percent
    scale_y_continuous(breaks = seq(0, 1, 0.20), limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 100, 20)) +
    # add text
    ylab("True/False Positive Rate") +
    labs(subtitle = "Boxplot by Method.",
         color = "Method") +
    facet_grid(~Method, scales='free_x', space='free', switch="x") +
    theme_bw(base_size = 24) +
    theme(strip.text.x = element_text(angle = 90, face="bold", colour="white"),
          strip.background = element_rect(fill="darkgray"),
          legend.title = element_text(size = 28),
          legend.text = element_text(size = 25)
    ) +
    labs(colour=guide_legend(title="Method"),
         caption = NULL,
         title = NULL,
         subtitle = NULL) +
    # xlab(str_glue("{descript}")) +
    guides(shape = "none")



geom_point(aes(shape = dataset), size = 4, alpha = 0.8) +
    facet_grid(~method_name, scales='free_x', space='free', switch="x") +
    theme_bw(base_size = 24) +
    geom_line(size = 1.9, alpha = 0.6)
### TPR Plots ----



load("analysis/robustness/Outputs/Cluster_Reshuffling/Iterator_Results_CR_seurat_pbmc_subset_top250_2021-12-10_14-55.RData")

load("analysis/robustness/Outputs/Cluster_Reshuffling/Iterator_Results_CR_seurat_pbmc_reshuffle_top100_2021-12-06_21-43.RData")
