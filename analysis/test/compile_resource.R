# Load required packages
require(tidyverse)
require(OmnipathR)
require(magrittr)
require(proxy)
require(viridis)
require(RCurl)
require(UpSetR)
require(liana)
require(shadowtext)
require(logger)
require(rlang)

# Source required functions
source("src/plot_utils.R")
source("analysis/resource_analysis/resource_descriptive.R")


## Resource descriptive analysis
# Obtain list with CCC Resources
ligrec <- compile_ligrec_descr()


# Run figure pipeline
descriptive_plots(ligrec)
saveRDS(.resource_env, "data/output/temp/resource_plots.RDS")

#
library(patchwork)
library(ggplotify)

#
resource_pdfs <- as.list(.resource_env)

resource_pdfs <- tibble(s_name = names(resource_pdfs),
                        plot = resource_pdfs %>% unname)





names(resource_pdfs)

c(
  # "upset", # 2x3 class = upset -> do manually?
  # these by classes and enrich
  "SignaLink",
  "SIGNOR",
  "NetPath",
  "CancerSEA",
  "MSigDB",
  "DisGeNet",
  "HPA_tissue_tissue",
  "HPA_tissue_organ",
  "HGNC",
  "OP-L"
  )

# c("jaccard", "shared") # 1x3
# Jaccard index plots
jaccs <- resource_pdfs %>%
    filter(str_detect(s_name, "jaccard")) %>%
    pluck("plot")

# shared
shared <- resource_pdfs %>%
    filter(str_detect(s_name, "shared")) %>%
    pluck("plot")

# ^ same for both
xd <- patchwork::wrap_plots(shared,
                            ncol=1,
                            nrow(3)) +
    plot_annotation(tag_levels = 'A',
                    tag_suffix = ')') &
    theme(plot.tag = element_text(face = 'bold',
                                  size = 32))
cairo_pdf(filename ="shared_comp.pdf",
          width = 16,
          height = 32
)
print(xd)
dev.off()


#
classes <- resource_pdfs %>%
    filter(str_detect(s_name, "classes"))


# by classes
classes <- resource_pdfs %>%
    filter(str_detect(s_name, "classes"))

#
classes



# Get LIGREC olap
ligrec %>%
    ligrec_decomplexify %T>%
    ligrec_overheats %>%
    ligrec_overlap %T>%
    uniq_per_res %T>%
    ligand_receptor_upset(upset_args = upset_args) %>%
    ligrec_classes_all %>%
    summarize_overlaps %T>%
    total_unique_bar



#
bench_env$`figures/descriptive/classes_genesets_CancerSEA.pdf`

s <- "figures/descriptive/classes_genesets_CancerSEA.pdf"
s <- gsub("(.*)\\/", "\\2", s) %>%
    gsub(".pdf*", "\\1", .)
s
