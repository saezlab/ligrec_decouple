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


### Mer/Seq Supp. Fig
### Seq/MerFISH
n_ranks <- c(50, 100,
             500, 1000,
             2500, 5000,
             10000)


seqfish_lr_coloc <- readRDS("data/output/spatial_out/fish/fish_lrcoloc.RDS") %>%
    filter(setting == "comp_n")


# 1) Initial Check
hist(seqfish_lr_coloc$estimate)
seqfish_lr_coloc %>%
    check_coloc()

# 2) FET boxplot
seqfish_lr_coloc %>%
    get_fet_boxplot_data(., n_ranks = n_ranks) %>%
    get_spatial_boxplot()
