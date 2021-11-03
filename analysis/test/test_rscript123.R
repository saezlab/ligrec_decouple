require(liana)
require(tibble)
require(purrr)
require(magrittr)

liana_path <- system.file(package = "liana")
testdata <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

liana_res <- liana_wrap(seurat_object, method = c('call_natmi',  'call_connectome',
                                                  'logfc', 'cellchat',
                                                  'call_sca', 'squidpy',
                                                  "cytotalk"))

saveRDS("data/output/test123.RDS")
