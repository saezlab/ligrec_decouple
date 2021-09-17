require(SingleCellExperiment)
require(Seurat)
require(liana)
require(tidyverse)
require(magrittr)




# make this thing into a function
# filter out any interactions that don't include the ADT-RNA intersect to save time. - tried gives different results :D
# perform it with 3 datasets + 3 intervals (0, 0.05, 0.1, 0.2)
# plot together
# analyse the proportions

# filter out control ADTs,
# look for alias genesymbols


# We Run Kendal Correlation, which supposedly uses tau-b scores (i.e. accounts for ties)
# Only 7 receptors match to OmniPath LRs, if I consider also ligands it would be a higher number,
# but very similar, but it would need to be separately, i.e. by source and target, so it becomes confusing
# and in my opinion redundant - ligand and receptors overlap also.


# Load matrix and convert to Seurat object ----
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

# Normalize ADT
seurat_object <- NormalizeData(seurat_object,
                               assay = "ADT",
                               normalization.method = "CLR")
# seurat_object <- ScaleData(seurat_object, assay = "ADT")


# save test seurat_objcet
saveRDS(seurat_object, "data/input/cmbc_seurat_test.RDS")



# RUN LIANA on cbmc test data ----
liana_res <- liana_wrap(seurat_object,
                        squidpy.params=list(cluster_key = "seurat_clusters",
                                            seed = as.integer(1)),
                        expr_prop = 0.2)
liana_res %<>% liana_aggregate()
saveRDS(liana_res, "data/output/test_citeseq_02prop.RDS")


# CiteSeq ----
## Read the above results
seurat_object <- readRDS("data/input/cmbc_seurat_test.RDS")
liana_res <- readRDS("data/output/test_citeseq_0prop.RDS")
op_resource <- select_resource("OmniPath")[[1]]

# convert to singlecell object
sce <- SingleCellExperiment::SingleCellExperiment(
    assays=list(counts = GetAssayData(seurat_object, assay = "ADT", slot = "counts"),
                data = GetAssayData(seurat_object, assay = "ADT", slot = "data")),
    colData=DataFrame(label=seurat_object@meta.data$seurat_clusters)
)



# check for alias----
# load the annotation database
library(org.Hs.eg.db)

# use sql to get alias table and gene_info table (contains the symbols)
# first open the database connection
dbCon <- org.Hs.eg_dbconn()
# write your SQL query
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# execute the query on the database
aliasSymbol <- DBI::dbGetQuery(dbCon, sqlQuery)

# set up your query genes
queryGeneNames <- rownames(sce)
# subset to get your results
map(queryGeneNames, function(name) aliasSymbol[which(aliasSymbol[,2] == stringr::str_to_upper(name)),5])



# Get ADT stats ----
# Get Summary per clust ----
test_summ <- scuttle::summarizeAssayByGroup(sce,
                                            ids = colLabels(sce),
                                            assay.type = "data")
test_summ@colData
means <- test_summ@assays@data$mean %>% # gene mean across cell types
    as_tibble(rownames = "entity_symbol")
means

props <- test_summ@assays@data$prop.detected %>% # gene prop
    as_tibble(rownames = "entity_symbol")
props

# Format Summary
receptor_syms <- c(op_resource$target_genesymbol)

adt_props_entity <- props %>%
    filter(entity_symbol %in% receptor_syms) %>%
    pivot_longer(-entity_symbol,
                 names_to = "cluster",
                 values_to = "adt_prop")

# means
adt_means_entity <- means %>%
    filter(entity_symbol %in% receptor_syms) %>%
    pivot_longer(-entity_symbol,
                 names_to = "cluster",
                 values_to = "adt_mean") %>%
    group_by(entity_symbol) %>%
    mutate(adt_scale = scale(adt_mean)) %>%
    unnest(adt_scale) %>%
    ungroup()

# join props to means
adt_means_entity %<>%
    left_join(adt_props_entity) %>%
    dplyr::rename("target" = cluster)


# Basic correlations between LRs and ADT means -----
liana_res2 <- liana_res %>%
    filter(receptor %in% adt_means_entity$entity_symbol) %>%
    select(source, ligand, target, receptor, ends_with("rank")) %>%
    pivot_longer(c(ligand,receptor),
                 names_to = "type",
                 values_to = "entity_symbol") %>%
    filter(type == "receptor") %>%
    distinct()

liana_adt <- liana_res2 %>%
    left_join(adt_means_entity) %>%
    select(ends_with("rank"), starts_with("adt")) %>%
    pivot_longer(-c(adt_scale, adt_mean, adt_prop))




# run corrs
set.seed(1)
adt_corr <- liana_adt %>%
    group_by(name) %>%
    nest() %>%
    mutate(mean_model = map(data, ~corr_fun(.x, var="adt_mean"))) %>%
    mutate(scale_model = map(data, ~corr_fun(.x, var="adt_scale"))) %>%
    mutate(prop_model = map(data, ~corr_fun(.x, var="adt_prop")))


# Check correlations
mean_corr <- adt_corr %>%
    select(name, mean_model) %>%
    unnest(cols = c(mean_model)) %>%
    mutate(mean_mlog10p = -log10(p.value)) %>%
    select(name, mean_pval = p.value, mean_estimate = estimate)

scale_corr <- adt_corr %>%
    select(name, scale_model) %>%
    unnest(cols = c(scale_model)) %>%
    select(name, scale_pval = p.value, scale_estimate = estimate)

prop_corr <- adt_corr %>%
    select(name, prop_model) %>%
    unnest(cols = c(prop_model)) %>%
    mutate(prop_mlog10p = -log10(p.value)) %>%
    select(name, prop_pval = p.value, prop_estimate = estimate)

corr_format <- mean_corr %>%
    # join
    left_join(scale_corr) %>%
    left_join(prop_corr) %>%
    dplyr::rename("method" = name) %>%
    pivot_longer(-method,
                 names_to = "corr_type") %>%
    separate(corr_type, into = c("metric", "stat")) %>%
    pivot_wider(names_from = stat,
                values_from = value) %>%
    mutate(significant = if_else(pval <= 0.05, TRUE, FALSE)) %>%
    select(-pval) %>%
    ungroup()

# plot
corr_format %>%
    # remove Mean/Median ranks
    filter(!(method %in% c("median_rank", "mean_rank"))) %>%
    # rename methods
    mutate(method = str_to_title(method)) %>%
    mutate(method = gsub("\\..*","",method)) %>%
    mutate(estimate = -1*estimate) %>% # reverse negative estimates
    # rename metrics
    mutate(metric = if_else(metric=="scale", "Cluster-specific Mean", metric)) %>%
    mutate(metric = if_else(metric=="mean", "Mean", metric)) %>%
    mutate(metric = if_else(metric=="prop", "Cell Proportion", metric)) %>%
    ggplot(aes(x = metric, y = estimate, colour = method, shape = significant)) +
    geom_point(size = 6, position = position_dodge(w = 0.05)) +
    xlab("Expression Type") +
    ylab("Kendal's tau Correlation Coefficients") +
    theme_minimal(base_size = 24) +
    labs(colour = "Method", shape="Significant")


# ROC according to Proportion thresholds ----
# hist(adt_means_entity[,-1] %>% as.matrix() %>% as.numeric(), breaks=100)
# normalize and cut anything below 0, i.e. assume that its noise

adt_props_entity
liana_prop <- liana_res2 %>%
    left_join(adt_props_entity, by = c("target" = "cluster", "entity_symbol"))





