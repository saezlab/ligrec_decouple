require(SingleCellExperiment)
require(Seurat)
require(liana)
require(tidyverse)
require(magrittr)

# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.1.0/5k_pbmc_protein_v3
mats <- Seurat::Read10X("data/input/10x_citeseq/5k_pbmc_protein_v3_filtered_feature_bc_matrix/filtered_feature_bc_matrix/")
rownames(mats$`Antibody Capture`) %<>%
    gsub("\\_.*","", .)



seurat_object <- SeuratObject::CreateSeuratObject(counts = mats$`Gene Expression`)
seurat_object[["ADT"]] <- CreateAssayObject(counts = mats$`Antibody Capture`)

# run def analysis
seurat_object %<>%
    FindVariableFeatures() %>%
    NormalizeData() %>%
    ScaleData() %>%
    RunPCA(verbose = TRUE) %>%
    FindNeighbors(reduction = "pca") %>%
    FindClusters(resolution = 0.4, verbose = TRUE)

# issue with squipy nums - need to fix
clust.anns <- str_glue("c{levels(seurat_object)}")
names(clust.anns) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, clust.anns)
seurat_object@meta.data <- seurat_object@meta.data %>%
    mutate(seurat_clusters = as.factor(str_glue("c{seurat_clusters}")))
Idents(seurat_object)

seurat_object

liana_res <- liana_wrap(seurat_object,
                        squidpy.params=list(cluster_key = "seurat_clusters"))
saveRDS(liana_res, "data/output/test_citeseq.RDS")


# Get ADT stats ----
# convert to singlecell object
sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts = GetAssayData(seurat_object, assay = "ADT")),
                                                  colData=DataFrame(label=seurat_object@meta.data$seurat_clusters))
op_resource <- select_resource("OmniPath")[[1]]


test_summ <- scuttle::summarizeAssayByGroup(sce,
                                            ids = colLabels(sce),
                                            assay.type = "counts")




test_summ@colData
means <- test_summ@assays@data$mean %>% # gene mean across cell types
    as_tibble(rownames = "entity_symbol")
means

props <- test_summ@assays@data$prop.detected %>% # gene prop
    as_tibble(rownames = "entity_symbol")
props

op_syms <- c(#op_resource$source_genesymbol, (only receptors)
             op_resource$target_genesymbol
             )

adt_props_entity <- props %>%
    filter(entity_symbol %in% op_syms) %>%
    pivot_longer(-entity_symbol,
                 names_to = "cluster",
                 values_to = "adt_prop")

# means
adt_means_entity <- means %>%
    filter(entity_symbol %in% op_syms)

# looks okay I guess, obviously need to be normalized before that
# scale and format
adt_scaled_entity <- adt_means_entity %>%
    pivot_longer(-entity_symbol,
                 names_to = "cluster",
                 values_to = "adt_mean") %>%
    group_by(entity_symbol) %>%
    mutate(adt_scale = scale(adt_mean)) %>%
    unnest(adt_scale) %>%
    ungroup()

# join props
adt_scaled_entity %<>%
    left_join(adt_props_entity) %>%
    dplyr::rename("target" = cluster)


# Basic corellations with LRs -----
liana_res <- liana_res %>% liana_aggregate()
liana_res2 <- liana_res %>%
    filter(receptor %in% adt_scaled_entity$entity_symbol) %>%
    select(source, ligand, target, receptor, ends_with("rank")) %>%
    pivot_longer(c(ligand,receptor),
                 names_to = "type",
                 values_to = "entity_symbol") %>%
    filter(type == "receptor") %>%
    distinct()

liana_adt <- liana_res2 %>%
    left_join(adt_scaled_entity) %>%
    select(ends_with("rank"), starts_with("adt")) %>%
    pivot_longer(-c(adt_scale, adt_mean, adt_prop))


#' Corr Function
#' @param df nested df
#' @param var name of the variable
corr_fun <- function(df, var){
    cor.test(df[[var]],
             df$value,
             method="kendall",
             exact = FALSE
    ) %>%
        broom::tidy()
}

# run corrs
adt_corr <- liana_adt %>%
    group_by(name) %>%
    nest() %>%
    mutate(mean_model = map(data, ~corr_fun(.x, var="adt_mean"))) %>%
    mutate(scale_model = map(data, ~corr_fun(.x, var="adt_scale"))) %>%
    mutate(prop_model = map(data, ~corr_fun(.x, var="adt_prop")))


# Check correlations
adt_corr %>%
    select(name, mean_model) %>%
    unnest(cols = c(mean_model))

adt_corr %>%
    select(name, scale_model) %>%
    unnest(cols = c(scale_model))

adt_corr %>%
    select(name, prop_model) %>%
    unnest(cols = c(prop_model))


mean_corr
scale_corr
