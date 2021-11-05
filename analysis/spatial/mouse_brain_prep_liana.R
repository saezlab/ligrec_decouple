# III) Run LIANA ----
# Load libs
require(tidyverse)
require(magrittr)
require(Seurat)
require(liana)

# brain analysis directory
brain_dir <- "data/input/spatial/brain_cortex"

cortex_sc <- readRDS(file.path(brain_dir, "allen_cortex_prep.rds"))
murine_resource <- readRDS("data/input/murine_omnipath.RDS")

liana_res <- liana_wrap(cortex_sc,
                        cellchat.params=list(organism="mouse"),
                        resource = "custom",
                        external_resource = murine_resource,
                        # this is passed only to squidpy, cellchat, cytotalk, and logfc
                        expr_prop=0.1,
                        cellchat.params = list(nboot=1000,
                                               expr_prop = 0),
                        call_natmi.params = list(
                            expr_file = str_glue("mouse_brain_em.csv"),
                            meta_file = str_glue("mouse_brain_metadata.csv"),
                            output_dir = str_glue("mouse_brain_results"),
                            reso_name = str_glue("mouse_brain_placeholder")
                        ),
                        assay = "SCT",
                        squidpy.params = list(cluster_key = "subclass"))


# squidpy sets gene names to upper (in the processing), revert this to title (i.e. murine)
# set all others to title just in case
liana_res %<>% map(function(res) res %>%
                       mutate_at(.vars = c("ligand", "receptor"), str_to_title))
saveRDS(liana_res, "data/input/spatial/brain_cortex/brain_liana_results.RDS")
