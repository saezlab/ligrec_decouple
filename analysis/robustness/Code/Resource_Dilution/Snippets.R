
# 1.
### Does the percentage of OmniPath genes in the results rise as OmniPath becomes more diluted? ----
# run after liana results after dilution are in env
{
gene_names_results <- unique(c(liana_results_OP$connectome$OmniPath_0$ligand, liana_results_OP$connectome$OmniPath_0$receptor))
gene_names_resource <- unique(c(resources_OP$connectome$OmniPath_0$source_genesymbol, resources_OP$connectome$OmniPath_0$target_genesymbol))
percentage_resource_in_results <- sum(gene_names_resource %in% gene_names_results)/length(gene_names_resource)
print("OmniPath_0")
print(percentage_resource_in_results)


for (i in names(dilution_props)) {
print(i)
gene_names_results <- unique(c(liana_results_OP$connectome[[i]]$ligand, liana_results_OP$connectome[[i]]$receptor))
gene_names_resource <- unique(c(resources_OP$connectome[[i]]$source_genesymbol, resources_OP$connectome[[i]]$target_genesymbol))

percentage_resource_in_results <- sum(gene_names_resource %in% gene_names_results)/length(gene_names_resource)
print(percentage_resource_in_results)

}

rm(i, gene_names_resource, gene_names_results, percentage_resource_in_results)
}


# 2. 
## does the number of interactions rise when lapplying over multiple OP resources regardless of dilution? ----
# run after dilutions are in env
{
test_list <- list(OmniPath_0 = resources_OP$connectome$OmniPath_0, 
                  OmniPath_0 = resources_OP$connectome$OmniPath_0, 
                  OmniPath_0 = resources_OP$connectome$OmniPath_0, 
                  OmniPath_0 = resources_OP$connectome$OmniPath_0)



test_result_list <- lapply(test_list, call_connectome, seurat_object = testdata)
lapply(test_result_list, all.equal, test_result_list$OmniPath_0)
## the results is the same all four times, and have the same number of rows
}

# 3. does the number of interactions in connectome rise when diluting only hits? ----
{
## if we filter OP_0 to be just hits and then dilute hits with more hits the number of rows, reflecting the number of interactions, shouldn't rise
# run after dilutions are in env
gene_names <- rownames(testdata@assays$RNA@data)
OmniPath_filter_0 <- resources_OP$connectome$OmniPath_0 %>%
                        filter(source_genesymbol %in% gene_names) %>%
                        filter(target_genesymbol %in% gene_names)
                
filter_connectome_0 <- call_connectome(seurat_object = testdata, op_resource = OmniPath_filter_0)
top_ranks_filter <- get_top_n_ranks(data_set =filter_connectome_0, method = "connectome", top_n = 200)

OmniPath_filter_60 <- dilute_Resource(resource = OmniPath_filter_0, top_rank_df = top_ranks_filter, dilution_prop = 0.6, data_set = testdata)
filter_connectome_60 <- call_connectome(seurat_object = testdata, op_resource = OmniPath_filter_60)
# the number of rows of filter_connectome_0 and filter_connectome_60 are the same


# as a test, we compare our connectome at filter_0 results with the connectome results of the standard
normal_connectome <- call_connectome(seurat_object = testdata, op_resource = resources_OP$connectome$OmniPath_0) %>%
  arrange_at(vars(everything()))
filter_connectome_0 <- call_connectome(seurat_object = testdata, op_resource = OmniPath_filter_0) %>%
  arrange_at(vars(everything()))

all.equal(normal_connectome, filter_connectome_0)
# the two results are identical





}

#4.
## Checking if a filtered (but undiluted) OP resource produces the same as the standard non-filtered one ----
## filtered means it only includes genes that are also present in the data you're going to use it with
## since gene interactions for genes not in the data should have no impact, this should make no difference.
## run when dilutions are in env
{
gene_names <- rownames(testdata@assays$RNA@data)
OmniPath_filter <- resources_OP$connectome$OmniPath_0 %>%
  filter(source_genesymbol %in% gene_names) %>%
  filter(target_genesymbol %in% gene_names)

filter_connectome <- call_connectome(op_resource = OmniPath_filter, seurat_object = testdata) %>%
                      arrange_at(vars(everything()))

all.equal(arrange_at(liana_results_OP$connectome$OmniPath_0, vars(everything())), filter_connectome)


filter_italk <- call_italk(op_resource = OmniPath_filter, seurat_object = testdata) %>%
  arrange_at(vars(everything()))

all.equal(arrange_at(liana_results_OP$italk$OmniPath_0, vars(everything())), filter_italk)

## works for italk and connectome when comparing to liana_wrap but not the others, but this is because liana_wrap and call_x don't work the same somehow
#### cellchat and sca are a bit different. I have to use the call function here but then the result is the same. See more in section below
filter_cellchat <- call_cellchat(op_resource = OmniPath_filter, seurat_object = testdata) %>%
  arrange_at(vars(everything()))
filter_cellchat_call <- call_cellchat(op_resource = select_resource(c("OmniPath"))[[1]], seurat_object = testdata) %>%
  arrange_at(vars(everything()))

all.equal(filter_cellchat_call, filter_cellchat)


filter_sca <- call_sca(op_resource = OmniPath_filter, seurat_object = testdata) %>%
  arrange_at(vars(everything()))
filter_sca_call <- call_sca(op_resource = select_resource(c("OmniPath"))[[1]], seurat_object = testdata) %>%
  arrange_at(vars(everything()))

all.equal(filter_sca_call, filter_sca)

}


#5.
#### Call_X and liana wrap don't produce the same results ----
{
#### even when same resource, same method, same data
## not the same output somehow
filter_sca_call <- call_sca(op_resource = select_resource(c("OmniPath"))[[1]], seurat_object = testdata)
filter_sca_wrap <- liana_wrap(seurat_object = testdata, method = c("sca"), resource = c("OmniPath"))[[1]]
## not the same output somehow
filter_cellchat_call <- call_cellchat(op_resource = select_resource(c("OmniPath"))[[1]], seurat_object = testdata)
filter_cellchat_wrap <- liana_wrap(seurat_object = testdata, method = c("cellchat"), resource = c("OmniPath"))[[1]]

}


# 6. Plot and save results for dilution, top 250 and top 1000 ranks -------------------------
{
load("C:/Users/plabu/OneDrive/Documents/GitHub/ligrec_robustness/Data/env_after_1000.RData")

top_rank_overlap_1000 <- top_rank_overlap_1000 %>%
  unnest(dilution_prop) %>%
  add_row(dilution_prop = 0, 
          connectome = 1, 
          cellchat = 1, 
          italk = 1, 
          sca = 1, 
          .before = 1) %>%
  "*"(100)



ggplot(data = top_rank_overlap_1000) + 
  geom_line(mapping = aes(dilution_prop, connectome, color =  "Connectome")) +
  geom_line(mapping = aes(dilution_prop, cellchat, color = "CellChat")) +
  geom_line(mapping = aes(dilution_prop, italk, color = "iTALK")) +
  geom_line(mapping = aes(dilution_prop, sca, color = "SCA")) +
  
  geom_point(mapping = aes(dilution_prop, connectome, color =  "Connectome")) +
  geom_point(mapping = aes(dilution_prop, cellchat, color = "CellChat")) +
  geom_point(mapping = aes(dilution_prop, italk, color = "iTALK")) +
  geom_point(mapping = aes(dilution_prop, sca, color = "SCA")) +
  
  ylim(0, 100) +
  
  ggtitle("Robustness of Method Predictions") +
  ylab("Overlap of Top Ranks [%]") +
  xlab("Dilution of Resource [%]") +
  labs(subtitle = "Generic dilution, top 1000 ranks",
       color = "Method")

ggsave("top_rank_overlap_1000_plot.png", 
       height = 5, width = 8, 
       path = "Outputs")

load("C:/Users/plabu/OneDrive/Documents/GitHub/ligrec_robustness/Data/env_after_250_high_res.RData")

top_rank_overlap_250 <- top_rank_overlap_250 %>%
  unnest(dilution_prop) %>%
  add_row(dilution_prop = 0, 
          connectome = 1, 
          cellchat = 1, 
          italk = 1, 
          sca = 1, 
          .before = 1) %>%
  "*"(100)



ggplot(data = top_rank_overlap_250) + 
  geom_line(mapping = aes(dilution_prop, connectome, color =  "Connectome")) +
  geom_line(mapping = aes(dilution_prop, cellchat, color = "CellChat")) +
  geom_line(mapping = aes(dilution_prop, italk, color = "iTALK")) +
  geom_line(mapping = aes(dilution_prop, sca, color = "SCA")) +
  
  geom_point(mapping = aes(dilution_prop, connectome, color =  "Connectome")) +
  geom_point(mapping = aes(dilution_prop, cellchat, color = "CellChat")) +
  geom_point(mapping = aes(dilution_prop, italk, color = "iTALK")) +
  geom_point(mapping = aes(dilution_prop, sca, color = "SCA")) +
  
  ylim(0, 100) +
  
  ggtitle("Robustness of Method Predictions") +
  ylab("Overlap of Top Ranks [%]") +
  xlab("Dilution of Resource [%]") +
  labs(subtitle = "Generic dilution, top 250 ranks",
       color = "Method")

ggsave("top_rank_overlap_250_plot.png", 
       height = 5, width = 8, 
       path = "Outputs")

}


## what is the unification of top ranks ----
{

unification_top_ranks <- top_ranks_OP_0$call_connectome[1:4] %>%
  add_row(top_ranks_OP_0$call_natmi[1:4]) %>%
  add_row(top_ranks_OP_0$call_italk[1:4]) %>%
  add_row(top_ranks_OP_0$call_sca[1:4]) %>%
  add_row(top_ranks_OP_0$cellchat[1:4])




}


## are tow result outputs identical when ordered? ----
{
all.equal(arrange_at(t, vars(everything())), 
          arrange_at(liana_results_OP$call_connectome$OmniPath_0, vars(everything()))
          )



}


## Investigating the relationship of OP source and target genesymbols ----
{
#get OP_resource and construct LR Pairs
op <- select_resource(c('OmniPath'))[["OmniPath"]] %>%
  select(source_genesymbol,
         target_genesymbol,
         is_directed,
         is_stimulation,
         consensus_stimulation,
         is_inhibition,
         consensus_inhibition,
         category_intercell_source,
         category_intercell_target,
         genesymbol_intercell_source,
         genesymbol_intercell_target,
         entity_type_intercell_target,
         sources,
         references,
         entity_type_intercell_source,
         entity_type_intercell_target) %>%
  mutate(isRandom = FALSE) %>%
  unite("LR_Pair", 
        c(source_genesymbol, target_genesymbol), 
        remove = FALSE, 
        sep = "_") %>%
  relocate("LR_Pair", .after = last_col())



# Are there duplicate Ligand-receptor pairs?
LR_Pair_list <- op$LR_Pair
length(unique(LR_Pair_list)) - length(LR_Pair_list) # No, every LR pair is unique

# Are there genes that are both source and target?
source_list <- op$source_genesymbol
target_list <- op$target_genesymbol

length(intersect(source_list, target_list)) # yes there are this many genes that occur as both source and target
length(unique(c(source_list, target_list))) # there are this many unique genesymbols in OP
length(intersect(source_list, target_list)) / length(unique(c(source_list, target_list))) # this proportion of genes in OP occur as both sources and targets


# Are any of these gene symbols source and target to themselves?
sum(source_list == target_list) # no, even though some genes appear in both lists they never appear twice in the same row, i.e. in relationship to themselves

}


# is thres = 1 actually a standard parameter of liana_wrap? ----
{
test_thing_no_thresh <- 
  liana_wrap(testdata, 
             method = c('cellchat'), 
             resource = c('OmniPath'), 
             expr_prop = 0,
             cellchat.params = list(nboot = cellchat_nperms, 
                                    expr_prop = 0),
             call_natmi.params = list(output_dir = natmi_output))

test_thing <- 
  liana_wrap(testdata, 
             method = c('cellchat'), 
             resource = c('OmniPath'), 
             expr_prop = 0,
             cellchat.params = list(nboot = cellchat_nperms, 
                                    expr_prop = 0,
                                    thresh = 1),
             call_natmi.params = list(output_dir = natmi_output))

}

# 11. NATMI produces more unique LR Pairs in the top ranking than other methods ----
# Other methods produce less unique LR Pairs and repeat them among more cluster combinations
# NATMI has less repeats
{

top_rank_edges <- 
  tibble("methods"             = script_params$methods_vector,
         "number_unique_edges" = 
           c(length(unique(top_ranks_OP$call_connectome$OmniPath_0$LR_Pair)),
             length(unique(top_ranks_OP$call_natmi$OmniPath_0$LR_Pair)),
             length(unique(top_ranks_OP$call_italk$OmniPath_0$LR_Pair)),
             length(unique(top_ranks_OP$call_sca$OmniPath_0$LR_Pair)),
             length(unique(top_ranks_OP$cellchat$OmniPath_0$LR_Pair))))

top_rank_edges <- top_rank_edges %>%
  arrange(desc(number_unique_edges))

print(top_rank_edges)




topology_tables <- list()
topology_plots  <- list()


for (method in top_rank_edges$methods) {

  # Investigate topology with plots
  topology_tables[[method]] <- top_ranks_OP[[method]]$OmniPath_0$LR_Pair %>%
    table()                                %>%
    as_tibble(.name_repair = "universal")  %>% 
    rename("LR_Pairs" = ".", "frequency" = "n")
  
  topology_tables[[method]]$LR_Pairs <- 
    reorder(topology_tables[[method]]$LR_Pairs, 
            topology_tables[[method]]$frequency)
  
  # often very large plot, adjuts plotting window accordingly
  topology_plots[[method]] <- 
    ggplot(data = topology_tables[[method]], 
           aes(x= LR_Pairs, y = frequency, width = 0.5)) +
    
    geom_bar(stat = "identity", fill = "#8ABAFF") +
    
    coord_flip() +
    
    ggtitle(str_glue("LR-Pair distribution in ", method)) +
    xlab("Unique Interactions") +
    ylab("Frequency among top-ranked interactions") +
    labs(subtitle = str_glue("Unique LR-Pairs : ",
                             nrow(topology_tables[[method]]),
                             " Highest frequency: ",
                             max(topology_tables[[method]]$frequency)))
  
  print(topology_plots[[method]])

}

rm(method)


# When you're done
rm(topology_plots, topology_tables, top_rank_edges)

}

# 12. Working env for updating dilute_Resource() ----
{
require(tidyverse)
require(Seurat)
require(liana)

library(lubridate)

# You may want to make sure the loaded function in the env are up to date
load("~/GitHub/ligrec_robustness/Data/Dilution_Upgrade_Env.RData")
master_seed <- 900

# This code to remove the old functions
rm(dilute_Resource, 
   get_top_n_ranks, 
   prop_isRandom, 
   rank_overlap, 
   extract_unconflicting_Genes, 
   preserve_Dilute, 
   random_Dilute)
}

# 13. Remove duplicates for resource_dilute [OBSOLETE] ----
{

# This code was part of the old heuristic version of random_Dilute(), it
# recursively removed duplicates and replaced them.
  
seed = 1

while(nrow(distinct(resource_dilute[, c("source_genesymbol", "target_genesymbol")])) < dilution_number) {
  
  
  print(str_glue("Removing Duplicates. Iteration: ", as.character(seed)))
  
  row_is_duplicated <- duplicated(resource_dilute[, c("source_genesymbol", 
                                                      "target_genesymbol")])
  
  source_gene_name_list <- gene_name_list
  
  set.seed(seed)
  resource_dilute[row_is_duplicated, ]$source_genesymbol <- 
    as.character(sample(source_gene_name_list,
                        size = sum(row_is_duplicated),
                        replace = TRUE))
  
  target_gene_name_list <- gene_name_list %>% 
    discard(~ .x %in% resource_dilute$source_genesymbol)
  
  
  set.seed(-seed)
  resource_dilute[row_is_duplicated, ]$target_genesymbol <- 
    as.character(sample(target_gene_name_list, 
                        size = sum(row_is_duplicated), 
                        replace = TRUE))
  
  
  seed = seed +1
  
  if(seed > 2000) {
    warning("MAXIMUM ITERATIONS REACHED. DUPLICATES REMOVAL FAILED.")
    break()
  }
}
  


}


# 14. Old code for liana dilutions call with separate top_rank_list
# Iterate over every method, lapply over every dilution proportion
for (method in methods_vector){
  
  dilutions_OP[[method]] <- 
    lapply(dilution_props, dilute_Resource, 
           resource          = resources_OP[[method]]$OmniPath_0, 
           top_rank_list     = top_ranks_OP[[method]]$OmniPath_0$LR_Pair, 
           preserve_topology = preserve_topology,
           data_set          = testdata,
           feature_type      = feature_type,
           verbose           = TRUE, 
           master_seed       = master_seed)
  
}



## 14. Check if all the testdata is the same ----

lapply(testdata, all.equal, testdata$Testdata_Seed_1)



## 15. Check if parts of the new results and unit test are the same

# unit test and new run agree on overlap
all_equal(robustness_old$complete_top_ranks_overlap[, c(1, 3, 4)], +
            robustness_standard$collated_top_ranks_overlap)

# unit test and new run don't agree on plot data structure but are visually 
# identical

# agreement liana results connectome 80% seed 1, 2
all_equal(
  robustness_old$dilution_robustness_results$liana_results_OP$call_connectome$OmniPath_80$OmniPath_80_Seed_1,
  robustness_standard$collated_robustness_results$liana_results_OP$call_connectome$OmniPath_80_Seed_1
)

all_equal(
  robustness_old$dilution_robustness_results$liana_results_OP$call_connectome$OmniPath_80$OmniPath_80_Seed_2,
  robustness_standard$collated_robustness_results$liana_results_OP$call_connectome$OmniPath_80_Seed_2
)


# same for natmi
all_equal(
  robustness_old$dilution_robustness_results$liana_results_OP$call_natmi$OmniPath_80$OmniPath_80_Seed_1,
  robustness_standard$collated_robustness_results$liana_results_OP$call_natmi$OmniPath_80_Seed_1
)

all_equal(
  robustness_old$dilution_robustness_results$liana_results_OP$call_natmi$OmniPath_80$OmniPath_80_Seed_2,
  robustness_standard$collated_robustness_results$liana_results_OP$call_natmi$OmniPath_80_Seed_2
)

# only the list names for testdata are different
all.equal(robustness_old$dilution_robustness_results$testdata, 
          robustness_standard$collated_robustness_results$testdata)

# resources are the same
all.equal(robustness_old$dilution_robustness_results$resources_OP$OmniPath_40$OmniPath_40_Seed_2, 
          robustness_standard$collated_robustness_results$resources_OP$OmniPath_40_Seed_2)

all.equal(robustness_old$dilution_robustness_results$resources_OP$OmniPath_80$OmniPath_80_Seed_1, 
          robustness_standard$collated_robustness_results$resources_OP$OmniPath_80_Seed_1)

all.equal(robustness_old$dilution_robustness_results$resources_OP$OmniPath_80$OmniPath_80_Seed_2, 
          robustness_standard$collated_robustness_results$resources_OP$OmniPath_80_Seed_2)




## 16. Setup testenvironment for editing the iterator




{
  # Load required Packages
  require(tidyverse)
  require(Seurat)
  require(liana)
  require(lubridate)
  
  
  # Define the functions needed to perform our analysis
  
  # Define the iterator wrapper function, which produces the end results
  source("Code/Resource_Dilution_Robustness/Robustness_Iterator.R")
  # Define Parameters for the following iterative robustness test
  source("Code/Resource_Dilution_Robustness/Iterator_Parameters.R")
  # Define functions needed for the iterator directly itself
  source("Code/Resource_Dilution_Robustness/Iterator_Functions.R")
  # Define functions that specifically reference the iterator parameters
  source("Code/Resource_Dilution_Robustness/Iterator_Parameter_Dependents.R")
  
  # Define functions for testing resource robustness
  source("Code/Resource_Dilution_Robustness/Resource_Dilution_Functions.R")
  # Define functions to dilute resources
  source("Code/Resource_Dilution_Robustness/Ranking_and_Misc_Functions.R")
  # Define functions to work with top_ranked CCIs and other miscellanea
  source("Code/Resource_Dilution_Robustness/Resource_Robustness_Functions.R")
  
  
  
  number_seeds      <- 2
  testdata_type     <- "liana_test"
  feature_type      <- "variable"
  preserve_topology <- FALSE
  dilution_props    <- c(seq(0.40, 1.00, 0.40))
  
  number_ranks <- list(
    "call_connectome" = 20,
    "squidpy"         = 20,
    "call_natmi"      = 20,
    "call_italk"      = 20,
    "call_sca"        = 20,
    "cellchat"        = 20
  )
  
  methods_vector <- c('call_connectome' ,
                      #'squidpy'         ,
                      #'call_natmi'      ,
                      'call_italk'      ,
                      'call_sca'        #,
                      #'cellchat'
  )
  
  sink_output     <- FALSE    
  liana_warnings  <- "divert" 
  
  save_results    <- TRUE
  trial_run       <- TRUE
  
  
  
  
  cellchat_nperms <- 10       
  
  bundled_outputs <- c(
    "liana_results_OP"  ,
    "resources_OP"      ,
    "top_ranks_OP"      ,
    "top_ranks_analysis",
    "runtime"          ,
    "testdata"
  )
  
  master_outputs <- c(
    "collated_top_ranks_overlap",
    "plot_box",
    "plot_line",
    "collated_robustness_results",
    "metadata"
  )                           
  
}


save.image("Outputs/Resource_Dilution/testenv_wrap_iterator.RData")


load("Outputs/Resource_Dilution/testenv_wrap_iterator.RData")






