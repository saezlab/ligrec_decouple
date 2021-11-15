require(tidyverse)
require(liana)
require(magrittr)

liana_path <- system.file(package = "liana")
testdata <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

liana_res <- liana_wrap(testdata, method = c('call_natmi',
                                             'call_connectome',
                                             'logfc', 'cellchat',
                                             'call_sca',
                                             'cellphonedb',
                                             "cytotalk"
                                             ),
                        call_natmi.params = list(reso_name = "placeholder123"),
                        permutation.params = list(parallelize=TRUE,
                                                  workers = 12),
                        resource = c("ICELLNET", "Default"))

saveRDS(liana_res, "data/output/test123.RDS")
