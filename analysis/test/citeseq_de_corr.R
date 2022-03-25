source("analysis/citeseq/citeseq_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")
source("analysis/comparison/comparison_utils.R")

require(liana)
require(tidyverse)
require(magrittr)
require(Seurat)
require(ggsignif)
require(yardstick)


# load prereqs
citeseq_dir <- "data/input/citeseq/"
arbitrary_thresh = 1.645 # one-tailed alpha = 0.05

# dataset <- "spleen_lymph_206"

ranks_of_interest <- c(50, 100, 250, 500, 1000, 2500, 5000)

datasets <- list.files(citeseq_dir)[c(1, 6, 7)]


de_corr_test <- map(datasets, function(dataset){
    message(dataset)

    # Load seurat
    seurat_object <- readRDS(file.path(citeseq_dir, dataset, str_glue("{dataset}_seurat.RDS")))
    liana_res <- readRDS(file.path(citeseq_dir, dataset, str_glue("{dataset}-liana_res-0.1.RDS")))


    # Set names to upper (needed for murine)
    rownames(seurat_object@assays$RNA@data) <- toupper(rownames(seurat_object@assays$RNA@data))
    rownames(seurat_object@assays$ADT@data) <- toupper(rownames(seurat_object@assays$ADT@data))

    # Obtain Barcodes by cluster
    barcode_by_cluster <- seurat_object@meta.data %>%
        rownames_to_column("barcodes") %>%
        dplyr::select(seurat_clusters, barcodes) %>%
        group_by(seurat_clusters)
    barcodes_key <- barcode_by_cluster %>%
        group_keys() %>%
        deframe()
    barcode_by_cluster %<>%
        group_split(.keep = FALSE) %>%
        setNames(barcodes_key) %>%
        map(function(x) deframe(x))


    # Generate aliases
    aliases <- get_adt_aliases(rownames(seurat_object@assays$ADT@data)) %>%
        mutate(across(everything(), ~toupper(.x)))

    # Keep only relevant symbols in each mat
    rna_mat <- seurat_object@assays$RNA@data[rownames(seurat_object@assays$RNA@data) %in% aliases$alias_symbol,]
    adt_mat <- seurat_object@assays$ADT@data[rownames(seurat_object@assays$ADT@data) %in% aliases$adt_symbol,]

    # Free RAM
    rm(seurat_object)
    gc()

    aliases


    corr_df <- map(barcode_by_cluster, function(barcodes){ # map across all clusters
        aliases %>%
            mutate(correlation = pmap(aliases, function(adt_symbol, # prot
                                                        alias_symbol # gene
            ){
                rna_vec <- rna_mat[rownames(rna_mat)==alias_symbol, colnames(rna_mat) %in% barcodes]
                adt_vec <- adt_mat[rownames(adt_mat)==adt_symbol, colnames(rna_mat) %in% barcodes]

                if(length(rna_vec)==0 | length(adt_vec)==0){
                    return(-999)
                }

                cor(rna_vec, adt_vec, method = "spearman")

            })) %>%
            unnest(correlation)
    }) %>%
        enframe(name = 'cluster', value = 'correlation') %>%
        unnest(correlation) %>%
        # remove NAs and missing vectors
        filter(correlation!=-999 & !is.na(correlation)) %>%
        # Z-transform by rna-prot pair across clusters
        group_by(adt_symbol, alias_symbol) %>%
        mutate(zcorr = scale(correlation)) %>%
        unnest(zcorr) %>%
        ungroup()

    corr_gt <- corr_df %>%
        dplyr::select(cluster, alias_symbol, zcorr)

    # aggregate liana
    liana_res <- liana_res %>%
        liana_aggregate(.score_mode = .score_comp, cap = 10000)

    liana_agg_gt <- liana_res %>%
        mutate(across(.cols=c("ligand", "receptor"), ~toupper(.x))) %>%
        dplyr::select(source, target, receptor, ends_with(".rank")) %>%
        dplyr::filter(receptor %in% corr_gt$alias_symbol) %>%
        mutate(across(c(source, target), ~str_replace_all(.x, " ", "."))) %>%
        mutate(across(ends_with("rank"), ~min_rank(.x))) %>% # convert to an integer rank
        pivot_longer(-c("source","target", "receptor"),
                     names_to = "method_name",
                     values_to = "predictor") %>%
        left_join(corr_gt, by=(c("target"="cluster",
                                 "receptor"="alias_symbol"))) %>%
        na.omit() %>%
        mutate(response = ifelse(abs(zcorr) >= arbitrary_thresh, 1, 0)) %>%
        mutate(response = factor(response, levels = c(1, 0)))

    ranked_results <- map(ranks_of_interest, function(n_rank){
        liana_agg_gt %>%
            group_by(method_name) %>%
            filter(predictor <= n_rank) %>%
            mutate(predictor = predictor * -1) %>%
            mutate(zcorr = abs(zcorr)) %>%
            mutate(correlation_test = cor(zcorr, predictor, method = "spearman")) %>%
            ungroup() %>%
            dplyr::select(method_name, correlation_test) %>%
            distinct() %>%
            mutate(n_rank = as.factor(n_rank)) %>%
            mutate(dataset = dataset)
    }) %>%
        bind_rows() %>%
        mutate(method_name = gsub("\\..*","", method_name)) %>%
        mutate(method_name = recode_methods(method_name))

    return(ranked_results)

})

de_corr_test %<>%
    bind_rows() %>%
    mutate(dataset = recode_datasets(dataset))

ggplot(de_corr_test,
       aes(x = n_rank, y = correlation_test,
           color = dataset, group=dataset)) +
    geom_point(aes(shape = dataset), size = 3, alpha = 0.8) +
    facet_grid(~method_name, scales='free_x', space='free', switch="x") +
    theme_bw(base_size = 21) +
    geom_line(size = 1.5, alpha = 0.6) +
    geom_hline(yintercept = 0, colour = "lightblue",
               linetype = 2, size = 1.7) +
    theme(strip.text.x = element_text(angle = 90, face="bold", colour="white"),
          axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5),
          strip.background = element_rect(fill="darkred"),
          legend.title = element_text(size = 23),
          legend.text = element_text(size = 21)
    ) +
    labs(colour=guide_legend(title="Dataset")) +
    ylab("Spearman Correlation") +
    xlab("Ranked Interactions Range") +
    guides(shape = "none")

