require(tidyverse)
require(magrittr)
require(liana)
require(RColorBrewer)
require(pheatmap)
require(proxy)
require(UpSetR)
require(grid)
require(ComplexHeatmap)
require(patchwork)
require(ggplotify)

source("analysis/comparison/comparison_utils.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

# Runs all things
# saves plot data
# patchwork
xxx <- comparison_pipe(input_filepath = "data/output/temp/liana_all_resources.RDS",
                       output_filepath = "test123",
                       resource = "OmniPath",
                       top_x = 0.05,
                       top_fun = "top_frac",
                       .score_specs = liana:::.score_specs,
                       cap_value_str = 1,
                       cap_value_freq = 1,
                       pval_thresh = 1,
                       sca_thresh = 0,
                       de_thresh = 1)

# Pre-defined resource list :) (max 8)
# e.g. > 2019 or something




# Save Supp Fig
output_filepath = "test123"
plot_path = "test123.pdf"
outpath <- as.character(str_glue("data/output/comparison_out/{output_filepath}/{plot_path}"))


cairo_pdf(outpath,
          width = 26,
          height = 32,
          family = 'DINPro')
(as.ggplot(xxx$jacc_heat) /
        (xxx$across_methods_jaccbox | xxx$across_resources_jaccbox)) +
    plot_layout(guides = 'collect', heights = c(3.5, 2.3)) +
    plot_annotation(tag_levels = 'A',
                    tag_suffix = ')') &
    theme(plot.tag = element_text(face = 'bold',
                                  size = 32))
dev.off()


cairo_pdf(outpath,
          width = 36,
          height = 50,
          family = 'DINPro')
(as.ggplot(xxx$freq_heat) /
        as.ggplot(xxx$strength_heat) /
        xxx$score_dist_plot) +
    plot_layout(guides = 'keep', heights = c(4, 4, 2.5)) +
    plot_annotation(tag_levels = 'A',
                    tag_suffix = ')') &
    theme(plot.tag = element_text(face = 'bold',
                                  size = 48),
          )
dev.off()


