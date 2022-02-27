# Load required packages
require(tidyverse)
require(OmnipathR)
require(magrittr)
require(proxy)
require(viridis)
require(RCurl)
require(liana)
require(shadowtext)
require(logger)
require(rlang)
library(patchwork)
library(ggplotify)
require(UpSetR)

# Source required functions
source("src/plot_utils.R")
source("analysis/resource_analysis/resource_descriptive.R")

## Resource descriptive analysis
# Obtain list with CCC Resources
ligrec <- compile_ligrec_descr()
saveRDS(ligrec, "data/input/ligrec.RDS")

# Run figure pipeline
ligrec <- readRDS("data/input/ligrec.RDS")
descript <- descriptive_plots(ligrec)
saveRDS(descript, "data/output/descript.RDS")

# Save Plots
descript <- readRDS("data/output/descript.RDS")

# Reproduce Compiled Plots
.resource_env <- descript
xtemp <- patchwork_resources()

# Assemble Main Figure 3 ----
path <- file.path("figures", "Main_Fig3.pdf")
pp <- patchwork::wrap_plots(list(
    as.ggplot(descript$size_overlap_combined +
                  theme(
                      text = element_text(size = 18),
                      axis.text.x = element_text(size=18, angle = 90,
                                                 vjust = 0.5, hjust=1),
                      axis.text.y = element_text(size=18)) +
                  xlab("Total Size")
              ),
    as.ggplot(descript$interactions_jaccard_heat)
    ),
    ncol=1,
    nrow(2)) +
    plot_layout(heights = c(1.3, 2)) +
    plot_annotation(tag_levels = list(c('A','B')), tag_suffix = ')') &
    theme(plot.tag = element_text(face = 'bold', size = 32))

cairo_pdf(filename = path,
          height = 17,
          width = 15)
print(pp)
dev.off()

# Assemble Main Fig. 4 ----
cust_theme <- theme(text = element_text(size = 13))

path <- file.path("figures", "Main_Fig4.pdf")

pp4 <- patchwork::wrap_plots(
    list(
        as.ggplot(descript$classes_perc_interactions_SignaLink_pathway + theme(text = element_text(size = 12),
                                                                                 legend.text = element_text(size=10))),
        as.ggplot(descript$enrich_heatmap_interactions_SignaLink_pathway + cust_theme),
        as.ggplot(descript$enrich_heatmap_interactions_CancerSEA + cust_theme),
        as.ggplot(descript$enrich_heatmap_interactions_HPA_tissue_organ + cust_theme)
        ),
    ncol=2, nrow(2)
    ) +
    plot_annotation(tag_levels = 'A', tag_suffix = ')') &
    theme(plot.tag = element_text(face = 'bold', size = 18))

cairo_pdf(filename = path,
          height = 7,
          width = 12.5)
print(pp4)
dev.off()



