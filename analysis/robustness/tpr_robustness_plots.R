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


robustness_plots <- imap(tpr_tobe_calculated, function(result_path, descript){
    if(stringr::str_detect(string = result_path, "Resource_Dilution")){
        analysis_focus = "resource"
    } else{
        analysis_focus = "cluster"
    }
    calculate_tpr(result_path, analysis_focus = analysis_focus) %>%
        plot_tpr(descript = descript)
})

require(patchwork)
path <- file.path( "figures", "SuppFig11_robustness.pdf")
pp <- patchwork::wrap_plots(
    robustness_plots,
    ncol=1,
    nrow(4)) +
    plot_annotation(tag_levels = 'A', tag_suffix = ')') &
    theme(plot.tag = element_text(face = 'bold', size = 32))
cairo_pdf(filename = path,
          height = 42,
          width = 18)
print(pp)
dev.off()
