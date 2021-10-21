require(tidyverse)
require(yardstick)
require(magrittr)
require(liana)
require(Seurat)
require(SeuratDisk)
library(ggsignif)
require(SingleCellExperiment)
require(decoupleR)
source("analysis/citeseq/citeseq_src.R")
source("src/eval_utils.R")


### Read Input ------

## read CytoSig
# models/signatures were obtained from https://github.com/data2intelligence/CytoSig/tree/master/CytoSig
cyto_signatures <- read.table("data/input/cytosig_signature_centroid.csv",
                              header = TRUE,
                              row.names = 1)

# Cytosig: add aliases (as in OP) and cytokine family members if appropriate
aliases_and_families <- list("CD40L" = "CD40LG",
                             "GSFC" = "CSF3",
                             "IFN1" = c("IFNA1", "IFNA2", "IFNA10",
                                        "IFNA7", "IFNA21", "IFNA5",
                                        "IFNA14", "IFNA17", "IFNA6",
                                        "IFNA4",  "IFNA16", "IFNA8" ),
                             "IFNL" = c("IFNL1", "IFNL2", "IFNL3", "IFNL4"),
                             "IL12" = c("IL12A", "IL12B"),
                             "IL36" = c("IL36A", "IL36B", "IL36G", "IL36RN"),
                             "MCSF" = "CSF1",
                             "TNFA" = c("TNF", "TNFA"),
                             "TRAIL" = c("TNFSF10"),
                             "TWEAK" = "TNFSF12"
)
alias_tib <- tibble(cytokine = names(aliases_and_families),
                    aliases = aliases_and_families) %>%
    unnest(aliases)


# Format Cytosig
cytosig_net <- cyto_signatures %>%
    as.data.frame() %>%
    rownames_to_column("target") %>%
    pivot_longer(-target,
                 names_to = "cytokine",
                 values_to = "weight") %>%
    # keep top 500 genes per cytokine
    group_by(cytokine) %>%
    slice_max(n=500, order_by = abs(weight)) %>%
    select(cytokine, target, weight) %>%
    mutate(mor = if_else(weight>0, 1, -1)) %>%
    mutate(weight = abs(weight)) %>%
    mutate(cytokine = if_else(cytokine=="A",
                              "Activin_A",
                              cytokine)) %>%
    ungroup()




#
# Calc Cytokine activities ----
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds")) %>%
    Seurat::NormalizeData()


require(decoupleR)
pseudobulk_counts <- Seurat::AverageExpression(seurat_object,
                                               assay = "RNA",
                                               slot = "data") %>%
    pluck("RNA")


# Run normalized mean
cytokine_enrichment <- run_wmean(
    pseudobulk_counts,
    cytosig_net,
    .source = "cytokine",
    .target = "target",
    .mor = "mor",
    .likelihood = "weight",
    times = 1000,
    seed = 42,
    sparse = TRUE,
    randomize_type = "cols_independently") %>%
    filter(statistic == "norm_wmean")

# Run normalized mean
cytokine_enrichment <- run_mlm(
    pseudobulk_counts,
    cytosig_net,
    .source = "cytokine",
    .target = "target",
    .mor = "mor",
    .likelihood = "weight",
    sparse = FALSE) %>%
    ungroup()


# Plot Barplot of Cytokine Enrichment
ggplot(cytokine_enrichment,
       aes(x = reorder(source, score), y = score)) +
    geom_bar(aes(fill = score), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred",
                         mid = "whitesmoke", midpoint = 0) +
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
          axis.text.x =
              element_text(angle = 45, hjust = 1, size =10, face= "bold"),
          axis.text.y = element_text(size =10, face= "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("Cytokines")

######
# Expand Cytosig results to aliases
cytokine_enrichment %<>%
    rename(cytokine = source) %>%
    left_join(alias_tib, by="cytokine") %>%
    mutate(aliases = if_else(is.na(aliases),
                             cytokine,
                             aliases)) %>%
    select(-cytokine, cytokine=aliases) %>%
    select(cytokine, everything())

# Check overlap
# Read OmniPath
op_resource <- select_resource("OmniPath")[[1]]

# Check which cytokines are not in OP
setdiff(unique(cytokine_enrichment$cytokine), unique(op_resource$source_genesymbol))

# Only keep CytoSig cytokines as ligands in OmniPath
op_resource %<>%
    filter(source_genesymbol %in% unique(cytokine_enrichment$cytokine))
op_resource

liana_res <- liana_wrap(seurat_object,
                        expr_prop = 0.1) %>%
    liana_aggregate()


liana_res_filt <- liana_res %>%
    filter(ligand %in% unique(cytokine_enrichment$cytokine))


#### Test with real data -------------------------------------------------------
assay.type <- "RNA" # maybe change to SCT?

seurat_object <- readRDS("data/input/spatial/Wu_etal_2021_BRCA/deconv/ER_celltype_minor/ER_celltype_minor_seurat.RDS")
Idents(seurat_object)
sce <- Seurat::as.SingleCellExperiment(seurat_object, assay = assay.type)
colLabels(sce) <- SeuratObject::Idents(seurat_object)

rm(seurat_object)
gc()

# Pseudobulk ----
pseudobulk <- scuttle::summarizeAssayByGroup(sce,
                                             ids = colLabels(sce),
                                             assay.type = "counts", # raw counts
                                             statistics = c("sum", "prop"))

# Subset metadata to only include the cluster and sample IDs to aggregate across
pseudobulk_expr <- pseudobulk@assays@data$sum %>%
    as_tibble(rownames = "gene") %>%
    pivot_longer(-gene, names_to = "celltype", values_to = "sum_count")
pseudobulk_prop <- pseudobulk@assays@data$prop.detected %>%
    as_tibble(rownames = "gene") %>%
    pivot_longer(-gene, names_to = "celltype", values_to = "prop")

pseudo <- pseudobulk_expr %>%
    left_join(pseudobulk_prop, by = c("gene", "celltype")) %>%
    # filter genes not expressed in at least 10% of cells per celltype
    # and only keep those with summed count of at least 10
    filter(prop > 0.1) %>%
    filter(sum_count >= 10) %>%
    select(-prop)

# looks good
pseudo %>%
    mutate(celltype = as.factor(celltype)) %>%
    ggplot(aes(x=celltype, y = log2(sum_count))) +
    geom_violin(trim=FALSE)

# Nest by Celltype, format, and normalize
pseudo %<>%
    group_by(celltype) %>%
    group_nest(.key = "counts") %>%
    # format and normalize
    mutate(logcounts = counts %>%
               map(function(c) c %>%
                       as.data.frame() %>%
                       column_to_rownames("gene") %>%
                       as.matrix() %>%
                       log2()))

# Enrichment
pseudo_cytosig <- pseudo %>%
    mutate(cytosig_res = logcounts %>%
               map(function(logc){
                   run_wmean(
                       logc,
                       cytosig_net,
                       .source = "cytokine",
                       .target = "target",
                       .mor = "mor",
                       .likelihood = "weight",
                       times = 1000,
                       seed = 1234,
                       sparse = TRUE,
                       randomize_type = "cols_independently") %>%
                       # keep only norm weighted mean
                       filter(statistic == "norm_wmean") %>%
                       # rename
                       select(cytokine=source,
                              NES=score,
                              p_value) %>%
                       # correct p
                       mutate(adj_pvalue = p.adjust(p_value))
                   }))
saveRDS(pseudo_cytosig, "data/output/er_brca_cytosig.RDS")

### Join LIANA to cytosig ====
pseudo_cytosig <- readRDS("data/output/er_brca_cytosig.RDS")
# cytosig format
cytosig_res <- pseudo_cytosig %>%
    select(celltype, cytosig_res) %>%
    unnest(cytosig_res) %>%
    left_join(alias_tib, by = "cytokine") %>%
    mutate(aliases = if_else(is.na(aliases),
                             cytokine,
                             aliases)) %>%
    select(-cytokine, cytokine=aliases) %>%
    select(cytokine, celltype, NES, adj_pvalue)

# read liana
liana_res <- readRDS("data/input/spatial/Wu_etal_2021_BRCA/deconv/ER_celltype_minor/ER_celltype_minor_liana_res.RDS") %>%
    liana_aggregate()


# join liana to cytosig
liana_res %<>%
    left_join(cytosig_res, by=c("source"="celltype"))

