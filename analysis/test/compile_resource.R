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
library(patchwork)
library(ggplotify)


# Source required functions
source("src/plot_utils.R")
source("analysis/resource_analysis/resource_descriptive.R")

#
upset_args <- list()


## Resource descriptive analysis
# Obtain list with CCC Resources
ligrec <- compile_ligrec_descr()


# Run figure pipeline
descriptive_plots(ligrec)
saveRDS(.resource_env, "data/output/temp/resource_plots.RDS")

# convert env to tibble
resource_pdfs <- tibble(s_name = names(as.list(.resource_env)),
                        plot = as.list(.resource_env) %>% unname)


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
          height = 32)
print(xd)
dev.off()

# Enrichments and Classes
dbs <- c("SignaLink",
         "SIGNOR",
         "NetPath",
         "CancerSEA",
         "MSigDB",
         "DisGeNet",
         "HPA_tissue_tissue",
         "HPA_tissue_organ",
         "HGNC",
         "OP-L")

# by class_enrich
sl_classes_enrich <- resource_pdfs %>%
    filter(str_detect(s_name, "classes_enrich"))

# by class_perc
sl_classes_perc <- resource_pdfs %>%
    filter(str_detect(s_name, "classes_perc"))



np_ce <- sl_classes_enrich %>%
    filter(str_detect(s_name, "NetPath")) %>%
    pluck("plot")


xd <- patchwork::wrap_plots(np_ce,
                            ncol=2,
                            nrow(3)) +
    plot_annotation(tag_levels = 'A',
                    tag_suffix = ')') &
    theme(plot.tag = element_text(face = 'bold',
                                  size = 12))
cairo_pdf(filename ="shared_comp.pdf",
          width = 12,
          height = 6)
print(xd)
dev.off()


# enrich
sl_classes_perc <- resource_pdfs %>%
    filter(str_detect(s_name, "enrich_heatmap"))



np_ce <- sl_classes_perc %>%
    filter(str_detect(s_name, "NetPath")) %>%
    pluck("plot")

patchwork_resources()


#' Helper Function to compile all resource plots
#'
#' @details makes use of `.resource_env` - a predefined environment to which I
#' use to save every resource plot. Then this function converts it into a list
#' and compiles all supp. figs.
patchwork_resources <- function(){
    # convert env to tibble
    resource_outs <- tibble(s_name = names(as.list(.resource_env)),
                            plot = as.list(.resource_env) %>% unname)

    # types of plots
    ptypes <- c("jaccard",
                "shared",
                "classes_enrich",
                "classes_perc",
                "enrich_heatmap")

    # external databases list
    dbs <- c("SignaLink",
             "SIGNOR",
             "NetPath",
             "CancerSEA",
             "MSigDB",
             "DisGeNet",
             "HPA_tissue_tissue",
             "HPA_tissue_organ",
             "HGNC",
             "OP-L")

    i <<- 0
    map(ptypes, function(plot_type){
        message(str_glue("Now compiling: {plot_type}"))

        # iterate for SuppFig Names
        i <<- i+1

        # filter by plot type
        resource_outs_filt <- resource_outs %>%
            filter(str_detect(s_name, plot_type))

        if(plot_type %in% c("jaccard", "shared")){

            path <- figure_path(
                'SuppFig_%s_%s.pdf',
                i, plot_type)

            pp <- patchwork::wrap_plots(resource_outs_filt %>% pluck("plot"),
                                        ncol=1,
                                        nrow(3)) +
                plot_annotation(tag_levels = 'A',
                                tag_suffix = ')') &
                theme(plot.tag = element_text(face = 'bold',
                                              size = 32))

            cairo_pdf(filename = path,
                      width = 16,
                      height = 32)
            print(pp)
            dev.off()


        } else if(plot_type %in% c("enrich_heatmap",
                                   "classes_enrich",
                                   "classes_perc")){

            map(dbs, function(db){
                path <- figure_path(
                    'SuppFig_%s_%s_%s.pdf',
                    i, plot_type, db)

                # additionally filter by db
                resource_outs_filt_plots <- resource_outs_filt %>%
                    filter(str_detect(s_name, db)) %>%
                    pluck("plot")

                # patchwork
                pp <- patchwork::wrap_plots(resource_outs_filt_plots,
                                            ncol=2,
                                            nrow(3)) +
                    plot_annotation(tag_levels = 'A',
                                    tag_suffix = ')') &
                    theme(plot.tag = element_text(face = 'bold',
                                                  size = 12))

                # to pdf
                cairo_pdf(filename = path,
                          width = 11,
                          height = 6)
                print(pp)
                dev.off()
            })
        }
    })
}





# by class_
# by class_perc
# by enrich
# + Supp. Fig Iterator


# Get LIGREC olap
ligrec_olap <- ligrec %>%
    ligrec_decomplexify %T>%
    ligrec_overheats %>%
    ligrec_overlap

#
ligrec_olap %T>%
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
