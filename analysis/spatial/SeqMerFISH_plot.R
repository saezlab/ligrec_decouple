############################## PLOT ############################################
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)

source("analysis/spatial/spatial_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

### Mer/Seq Supp. Fig
### Seq/MerFISH
n_ranks <- c(50, 100,
             250, 500,
             1000, 2500,
             5000)


seqfish_lr_coloc <- readRDS("data/output/spatial_out/fish/fish_lrcoloc.RDS") %>%
    filter(setting == "comp_n") %>%
    na.omit()

# 1) Initial Check
hist(seqfish_lr_coloc$estimate)
seqfish_lr_coloc %>%
    check_coloc()

# 2) FET boxplot
p <- seqfish_lr_coloc %>%
    get_fet_boxplot_data(., n_ranks = n_ranks) %>%
    get_eval_boxplot(eval_type = "space")


# Print Supp Fig
path <- file.path("figures",
                  "SuppFigure19_Fish_Spatial.RDS")
cairo_pdf(path,
          height = 12,
          width = 18,
          family = 'DINPro')
print(p)
dev.off()






