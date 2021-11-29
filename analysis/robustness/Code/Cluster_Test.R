require(liana)
require(Seurat)
require(tidyverse)
require(lubridate)

liana_path <- system.file(package = "liana")
testdata <-
  readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
# testdata <-  readRDS("Data/pbmc3k_final.rds")

liana_results <- liana_wrap(testdata,
                            method = c('call_connectome',
                                        'call_natmi',
                                        'call_italk',
                                        'call_sca',
                                        'cellchat',
                                        'squidpy'),
                            resource = c("OmniPath"))


test_plot <- ggplot(liana_results$cellchat,
                    aes(x = pval)) +
  geom_histogram(bins = 20)

ggsave(plot = test_plot, 
       "cluster_test.png", 
       height = 8.00,
       width = 8.00,
       path = "Outputs")

print("Hello world!")

print(Sys.time())

warning("uh oh...")

save(liana_results, file = "Outputs/cluster_test.RData")
