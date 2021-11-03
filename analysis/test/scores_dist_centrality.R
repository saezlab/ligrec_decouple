# LOAD New Results from different Method-Resource Combinations
liana_all <- readRDS("data/output/liana_all_resources.RDS")


liana_res <- liana_wrap(seurat_object, method = c('call_natmi',  'call_connectome',
                                                  'logfc', 'cellchat',
                                                  'call_sca', 'squidpy',
                                                  "cytotalk"))
