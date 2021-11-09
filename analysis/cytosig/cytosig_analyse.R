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

# path_tibble (relative paths to the relevant objects)
path_tibble <- tibble(dataset = c("ER",
                                  "HER2",
                                  "TNBC"),
                      seurat_path = c("data/input/spatial/Wu_etal_2021_BRCA/deconv/ER_celltype_minor/ER_celltype_minor_seurat.RDS",
                                      "data/input/spatial/Wu_etal_2021_BRCA/deconv/HER2_celltype_minor/HER2_celltype_minor_seurat.RDS",
                                      "data/input/spatial/Wu_etal_2021_BRCA/deconv/TNBC_celltype_minor/TNBC_celltype_minor_seurat.RDS"
                                      ),
                      liana_path = c(file.path("data/output/comparison_out/", str_glue("BRCA_TNBC_liana_OmniPath.RDS")),
                                     file.path("data/output/comparison_out/", str_glue("BRCA_HER2_liana_OmniPath.RDS")),
                                     file.path("data/output/comparison_out/", str_glue("BRCA_TNBC_liana_OmniPath.RDS"))
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
                                                     generate = TRUE) #!
                        gc()
                        return(cyto_res)
    }))
saveRDS(cytosig_eval, "data/output/cytosig_out/brca_cytosig_res.RDS")

# Read results
cytosig_eval <- readRDS("data/output/cytosig_out/brca_cytosig_res.RDS") %>%
    select(dataset, cytosig_res) %>%
    unnest(cytosig_res)

aucs <- cytosig_eval %>%
    select(-c(cyto_liana, corr)) %>%
    unnest(prc) %>%
    select(dataset, method_name, roc, prc_auc = auc) %>%
    distinct() %>%
    unnest(roc) %>%
    select(dataset, method_name, roc_auc = auc, prc_auc) %>%
    distinct() %>%
    group_by(method_name) %>%
    mutate(roc_mean = mean(roc_auc),
           prc_mean = mean(prc_auc)) %>%
    ungroup() %>%
    mutate(method_name = gsub("\\..*","", method_name)) %>%
    mutate(method_name = recode_methods(method_name))


roc_min <- ifelse(min(aucs$roc_auc) > 0.5, 0.5, min(aucs$roc_auc))
prc_min <- ifelse(min(aucs$prc_auc) > 0.5, 0.5, min(aucs$prc_auc))

# min_lim <- floor(min(c(aucs$roc, aucs$prc)) * 100)/100
# max_lim <- ceiling(max(c(aucs$roc, aucs$prc)) * 100)/100

ggplot(aucs,
       aes(x=roc_mean,
           y=prc_mean,
           color=method_name)) +
    geom_point(shape = 9, size = 12, alpha=1) +
    geom_point(aes(x = roc_auc,
                   prc_auc,
                   shape=dataset),
               size = 6,
               alpha = 0.3) +
    theme(text = element_text(size=16)) +
    xlab('AUROC') +
    ylab('AUPRC') +
    xlim(roc_min, 0.7) +
    ylim(prc_min, 0.8) +
    geom_hline(yintercept = 0.5, colour = "pink",
               linetype = 2, size = 1.2) +
    geom_vline(xintercept = 0.5, colour = "pink",
               linetype = 2, size = 1.2) +
    theme_bw(base_size = 30) +
    guides(shape=guide_legend(title="Dataset"),
           color=guide_legend(title="Method"))

