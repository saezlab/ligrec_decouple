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
cytosig_eval_wrap(.eval = "independent",
                  score_mode = "specs",
                  generate = TRUE) # only generate cytosig enrichments once

# ## B) Missing imputed as max
cytosig_eval_wrap(.eval = "max",
                  score_mode = "specs",
                  generate = FALSE)

## C) Intersect
cytosig_eval_wrap(.eval = "intersect",
                  score_mode = "specs",
                  generate = FALSE)


# II) Mixed/Comp Scores ----
## A) Independent (missing imputed as NA)
cytosig_eval_wrap(.eval = "independent",
                  score_mode = "mixed",
                  generate = FALSE)

## B) Missing imputed as max
cytosig_eval_wrap(.eval = "max",
                  score_mode = "mixed",
                  generate = FALSE)

## C) Intersect
cytosig_eval_wrap(.eval = "intersect",
                  score_mode = "mixed",
                  generate = FALSE)


# III) Housekeeping Scores
## A) Independent (missing imputed as NA)
cytosig_eval_wrap(.eval = "independent",
                  score_mode = "house",
                  generate = FALSE)

# ## B) Missing imputed as max
cytosig_eval_wrap(.eval = "max",
                  score_mode = "house",
                  generate = FALSE)

## C) Intersect
cytosig_eval_wrap(.eval = "intersect",
                  score_mode = "house",
                  generate = FALSE)
