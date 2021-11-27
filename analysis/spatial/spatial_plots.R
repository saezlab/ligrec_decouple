############################## PLOT ############################################
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)

source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

# Boxplot (Mixed)
cairo_pdf(file.path("data/output/temp/",
                    str_glue("spatial_fets_mixed.RDS")),
          height = 16,
          width = 24,
          family = 'DINPro')
print(get_spatial_boxplot("data/output/spatial_out/all_fets_mixed.RDS"))
dev.off()


# Boxplot (Mixed)
cairo_pdf(file.path("data/output/temp/",
                    str_glue("spatial_fets_specs.RDS")),
          height = 16,
          width = 24,
          family = 'DINPro')
print(get_spatial_boxplot("data/output/spatial_out/all_fets_specs.RDS"))
dev.off()

