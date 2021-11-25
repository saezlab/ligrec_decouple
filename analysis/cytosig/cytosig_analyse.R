# load required libraries
require(magrittr)
require(tidyverse)
require(liana)
require(Seurat)
require(yardstick)
require(SingleCellExperiment)
require(decoupleR)

source("analysis/cytosig/cytosig_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

# Check cytosig_eval_wrap for details!

# I) Specs ----
## A) Independent (missing imputed as NA)
cytosig_eval_wrap(eval = "independent",
                  score_mode = "specs",
                  generate = TRUE)

## B) Missing imputed as max
cytosig_eval_wrap(eval = "max",
                  score_mode = "specs",
                  generate = FALSE)

## C) Intersect
cytosig_eval_wrap(eval = "intersect",
                  score_mode = "specs",
                  generate = FALSE)


# II) Mixed/Comp Scores ----
cytosig_eval_wrap(eval = "independent",
                  score_mode = "mixed",
                  generate = FALSE)

## B) Missing imputed as max
cytosig_eval_wrap(eval = "max",
                  score_mode = "mixed",
                  generate = FALSE)

## C) Intersect
cytosig_eval_wrap(eval = "intersect",
                  score_mode = "mixed",
                  generate = FALSE)


# III) Housekeeping Scores
cytosig_eval_wrap(eval = "independent",
                  score_mode = "house",
                  generate = FALSE)

## B) Missing imputed as max
cytosig_eval_wrap(eval = "max",
                  score_mode = "house",
                  generate = FALSE)

## C) Intersect
cytosig_eval_wrap(eval = "intersect",
                  score_mode = "house",
                  generate = FALSE)




# path_tibble (relative paths to the relevant objects)
path_tibble <- tibble(dataset = c("ER",
                                  "HER2",
                                  "TNBC"),
                      # BRCA seurat objects used throughout the manuscript
                      seurat_path = c("data/input/spatial/Wu_etal_2021_BRCA/deconv/ER_celltype_minor/ER_celltype_minor_seurat.RDS",
                                      "data/input/spatial/Wu_etal_2021_BRCA/deconv/HER2_celltype_minor/HER2_celltype_minor_seurat.RDS",
                                      "data/input/spatial/Wu_etal_2021_BRCA/deconv/TNBC_celltype_minor/TNBC_celltype_minor_seurat.RDS"
                      ),
                      # liana results generated from extract_evals.sh
                      liana_path = c(file.path("data/output/brca_extracts", str_glue("ER_{eval}_{score_mode}_liana_res.RDS")),
                                     file.path("data/output/brca_extracts/", str_glue("HER2_{eval}_{score_mode}_liana_res.RDS")),
                                     file.path("data/output/brca_extracts", str_glue("TNBC_{eval}_{score_mode}_liana_res.RDS"))
                      ))

# get cytosig
cytosig_net <- load_cytosig()

# Run Cytosig Evaluation
cytosig_eval <- path_tibble %>%
    mutate(cytosig_res =
               pmap(path_tibble,
                    function(dataset, seurat_path, liana_path){
                        # Load Seurat Object
                        message(str_glue("Now running: {dataset}"))

                        seurat_object <- readRDS(seurat_path)
                        print(seurat_object)

                        # read liana
                        liana_res <- readRDS(liana_path)

                        # Run Cytosig Evaluation
                        cyto_res <- run_cytosig_eval(seurat_object = seurat_object,
                                                     liana_res = liana_res,
                                                     cytosig_net = cytosig_net,
                                                     z_scale = FALSE,
                                                     expr_prop = 0.1,
                                                     assay = "RNA",
                                                     sum_count_thresh = 5,
                                                     NES_thresh = 1.645,
                                                     subtype = dataset,
                                                     generate = generate) #!
                        gc()
                        return(cyto_res)
                    }))
saveRDS(cytosig_eval, "data/output/cytosig_out/cytosig_res_{eval}_{score_mode}.RDS")

