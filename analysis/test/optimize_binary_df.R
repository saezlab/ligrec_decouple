## Script to generate the comparison data presented in the Manuscript
library(tidyverse)
library(magrittr)
library(liana)
library(RColorBrewer)
library(pheatmap)
library(proxy)
library(UpSetR)
library(grid)
library(ComplexHeatmap)
library(patchwork)
library(ggplotify)

source("analysis/comparison/comparison_utils.R")
source("src/eval_utils.R")
source("src/plot_utils.R")

# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata",
                      "input", "testdata.rds"))

liana_res <- liana_wrap(seurat_object,
                        method = c("call_connectome", "connectome"),
                        resource = c("OmniPath", "ICELLNET"))
liana_agg <- liana_res %>% liana_aggregate(resource = "ICELLNET")

saveRDS(liana_res, "data/output/temp/test_binary_df.RDS")

# Ranked Scores according to a set of criteria (here by specificy whenever available)
liana_all_spec <- get_spec_list("data/output/temp/test_binary_df.RDS",
                                .score_spec = liana:::.score_specs)
liana_all_spec$connectome@method_results$OmniPath %<>%
    mutate(p_val_adj.rec = 0,
           p_val_adj.lig = 0)
liana_all_spec$connectome@method_results$ICELLNET %<>%
    mutate(p_val_adj.rec = 0,
           p_val_adj.lig = 0)


top_lists <- get_top_hits(liana_all_spec,
                          n_ints = 500,
                          top_fun = "top_n",
                          de_thresh = 1
                          )

jaccard_per_mr <- simdist_resmet(top_lists$top_500,
                                 sim_dist = "simil",
                                 method = "Jaccard")
jaccard_per_mr$reso
# ^ Works as expected


# Re-work get_binary
sig_list <- top_lists$top_500
binary_df <- get_binary_df(sig_list)


# get method and resource names combined
lnames <- map(names(sig_list), function(m_name){
    map(names(sig_list[[m_name]]), function(r_name){
        str_glue("{m_name}⊎{r_name}")
    })
}) %>%
    unlist()

# get binarized significant hits list (1 for sig per method, 0 if absent)
named_list <- sig_list %>%
    purrr::flatten() %>%
    setNames(lnames)
named_list %>% prepForUpset

named_list_binarized <- map(names(named_list), function(l_name){

    if(nrow(named_list[[l_name]])==0){
        message(str_glue("{l_name} has 0 hits"))
        return(tibble(interaction=character()))
    }

    named_list[[l_name]] %>%
        select(source, target, ligand, receptor) %>%
        unite("interaction", source, target,
              ligand, receptor, sep="_") %>%
        mutate(!!l_name := 1) %>%
        distinct() %>%
        data.table::setDT() %>%
        data.table::setkey(., interaction)
})


merge_test <- function(){
    named_list_binarized %>%
        reduce(., merge, by = "interaction") %>%
        mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
        mutate_at(vars(2:ncol(.)), ~ replace(., . != 0, 1))
}

full_join_test <- function(){
    named_list_binarized %>%
        reduce(., full_join, by = "interaction") %>%
        mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
        mutate_at(vars(2:ncol(.)), ~ replace(., . != 0, 1)) %>%
        as.data.frame()
}


MergedDT <- Reduce(function(...) merge(..., all = TRUE), named_list_binarized) %>%
    as.data.frame() %>%
    mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
    mutate_at(vars(2:ncol(.)), ~ replace(., . != 0, 1)) %>%
    column_to_rownames("interaction")


named_list_binarized %>%
    reduce(., left_join, by = "interaction")

named_list_binarized %>%
    reduce(., right_join, by = "interaction")


peakRAM(MergedDT = Reduce(function(...) merge(..., all = TRUE), named_list_binarized))
peakRAM(get_binary_df(sig_list))


all_equal(xd, MergedDT)


xd <- get_binary_df(sig_list)
xd2 <- get_binary_df(sig_list)

