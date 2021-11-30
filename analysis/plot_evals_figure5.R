# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(patchwork)

source("analysis/cytosig/cytosig_src.R")
source("src/plot_utils.R")
source("analysis/spatial/spatial_src.R")

# Mixed/Comp Method specifics
space_mixed <- get_spatial_bigbox("data/output/spatial_out/all_fets_mixed.RDS")
cytosig_mixed <- plot_cytosig_aucs(.eval = "independent",
                                   score_mode = "mixed")

path <- file.path("figures",
                  "Figure5_Evals_Composite.RDS")
cairo_pdf(path,
          height = 26,
          width = 22,
          family = 'DINPro')
print((space_mixed /
        cytosig_mixed) +
    plot_layout(guides = 'keep', heights = c(1, 1)) +
    plot_annotation(tag_levels = 'A',
                    tag_suffix = ')') &
    theme(plot.tag = element_text(face = 'bold',
                                  size = 40)))
dev.off()


# Specificities Method specifics
space_specs <- get_spatial_bigbox("data/output/spatial_out/all_fets_specs.RDS")
cytosig_specs <- plot_cytosig_aucs(.eval = "independent", score_mode = "specs")

path <- file.path("figures",
                  "SuppFigure19_Evals_Specific.RDS")
cairo_pdf(path,
          height = 26,
          width = 22,
          family = 'DINPro')
print((space_specs /
           cytosig_specs) +
          plot_layout(guides = 'keep', heights = c(1, 1)) +
          plot_annotation(tag_levels = 'A',
                          tag_suffix = ')') &
          theme(plot.tag = element_text(face = 'bold',
                                        size = 40)))
dev.off()
