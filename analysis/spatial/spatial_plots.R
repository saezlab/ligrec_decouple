############################## PLOT ############################################
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)

source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")




type_lr_coloc <- readRDS("data/output/spatial_out/all_fets_specs.RDS")
# Boxplot
cairo_pdf(file.path("data/output/temp/",
                    str_glue("spatial_fets_specs.RDS")),
          height = 24,
          width = 32,
          family = 'DINPro')
print(
    ggplot(type_lr_coloc,
           aes(x = n_rank,
               y = odds_ratio,
               color = dataset_type,
               fill = dataset_type)) +
        geom_boxplot(alpha = 0.2,
                     outlier.size = 1.5,
                     width = 0.8)  +
        # geom_point(aes(shape = dataset), size = 2) +
        scale_shape_manual(values = rep(1:12, len =  length(unique(type_lr_coloc$dataset)))) +
        facet_grid(~method_name, scales='free_x', space='free', switch="x") +
        theme_bw(base_size = 24) +
        geom_hline(yintercept = 1, colour = "black",
                   linetype = 2, size = 1.3) +
        theme(strip.text.x = element_text(angle = 90),
              axis.text.x = element_text(angle = 45, hjust=1)
        ) +
        # guides(color = "none") +
        labs(shape = guide_legend(title="Visium slide"),
             color = guide_legend(title="Dataset type")) +
        guides(fill = "none") +
        ylab("Odds Ratio") +
        xlab("#Ranks Considered")
)
dev.off()



