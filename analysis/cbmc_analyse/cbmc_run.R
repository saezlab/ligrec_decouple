# Run on local
require(liana)
require(tidyverse)
require(magrittr)
require(Seurat)

# Load and and format Data
# SeuratData::InstallData("cbmc")

cbmcdata <- SeuratData::LoadData("cbmc")
cbmcdata@meta.data %<>%
    mutate(celltype = as.factor(rna_annotations)) %>%
    # remove + as it breaks Squidpy
    mutate(rna_annotations = str_replace_all(rna_annotations, "[+]", "")) %>%
    mutate(rna_annotations = str_replace(rna_annotations, "/", " "))  %>%
    mutate(rna_annotations = if_else(rna_annotations=="B", "B cell", rna_annotations)) %>%
    # filter multiplets and mouse droplets
    filter(!(rna_annotations %in% c("Mouse", "Multiplets", "T Mono doublets")))
cbmcdata <- subset(cbmcdata, cells = rownames(cbmcdata@meta.data))
cbmcdata <- Seurat::SetIdent(cbmcdata, value = cbmcdata@meta.data$rna_annotations)
seurat_object <- cbmcdata %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures()
# name to seurat_cluster for consistency with other cmbc seurat objects
seurat_object@meta.data %<>%
    mutate(seurat_clusters = as.factor(rna_annotations))
Idents(seurat_object) <- seurat_object@meta.data$seurat_clusters
saveRDS(seurat_object, "data/input/citeseq/cmbcs/cbmc_seurat.RDS")

seurat_object <- readRDS("data/input/cbmc_seurat.RDS")




# Call liana
cmbc_results <- liana_wrap(seurat_object = seurat_object,
                           resource = "OmniPath",
                           method = c("natmi", "connectome", "logfc",
                                      "sca"),
                           expr_prop = 0,
                           squidpy.param = list(cluster_key = "celltype"))
cmbc_aggregate <- cmbc_results %>%
    liana_aggregate()


# Get ADT assay
adt_assay <- GetAssayData(seurat_object, slot = "data", assay = "ADT")
GetAssayData(seurat_object, assay = "ADT")


# check ligs and receptors that actually match the adts
ligs <- unique(cmbc_aggregate$ligand)
recpts <- unique(cmbc_aggregate$receptor)
ligs_present <- ligs[ligs %in% rownames(adt_assay)]
recpts_present <- recpts[recpts %in% rownames(adt_assay)]


# extract only relevant lrs
cmbc_relevant_lrs <- cmbc_aggregate %>%
    filter(receptor %in% c(ligs_present, recpts_present)) %>%
    pivot_longer(cols=c(ligand, receptor),
                 names_to = "entity",
                 values_to = "entity_symbol") %>%
    filter(entity_symbol %in% c(ligs_present, recpts_present)) %>%
    select(source, target, entity_symbol, ends_with("rank"))


# extract relevant ADT assay expression
sce <- Seurat::as.SingleCellExperiment(seurat_object,  assay = "ADT")
SingleCellExperiment::colLabels(sce) <- SeuratObject::Idents(seurat_object)
# sce@assays@data$scaledata <- seurat_object@assays$RNA@scale.data

mean_prop <-
    scuttle::summarizeAssayByGroup(sce,
                                   ids = SingleCellExperiment::colLabels(sce),
                                   assay.type = "counts",
                                   statistics = c("mean", "prop.detected"))


mean_prop


# seurat_object <- NormalizeData(seurat_object,
#                                normalization.method = "CLR",
#                                assay= "ADT")
avg_expr <- Seurat::AverageExpression(seurat_object) %>%
    pluck("ADT")

adt_values <- avg_expr %>%
    as.tibble(rownames = "entity_symbol") %>%
    filter(entity_symbol %in% c(ligs_present, recpts_present)) %>%
    pivot_longer(-entity_symbol,
                 values_to = "adt_expression",
                 names_to = "celltype")




# join lrs and adt expression
adt_lrs <- cmbc_relevant_lrs %>%
    left_join(adt_values %>%
                  rename(source = celltype) %>%
                  rename(adt_source = adt_expression),
              by=c("entity_symbol", "source")) %>%
    left_join(adt_values %>%
                  rename(target = celltype) %>%
                  rename(adt_target = adt_expression),
              by=c("entity_symbol", "target"))

# correlation to threshold expressing these things
# specificty - ADTs are not expressed
plot(density(adt_lrs %>% filter(entity_symbol == "CD4") %>% pull("adt_target") %>% unique()))


# cmbc_cc <- liana_wrap(
#     method = c("connectome"),
#     seurat_object = seurat_object,
#     resource = c("OmniPath", "CellPhoneDB"),
#     natmi.params = list(
#         expr_file = "em.csv",
#         meta_file = "metadata.csv",
#         output_dir = "NATMI_cbmc",
#         assay = "RNA",
#         num_cor = 4,
#         .format = TRUE,
#         .write_data = FALSE,
#         .seed = 1004,
#         .natmi_path = "~/Repos/LIANA/NATMI"),
#     squidpy.params = list(
#         cluster_key="rna_annotations",
#         n_perms=10000,
#         threshold=0.1,
#         seed=as.integer(1004)
#     ),
#     cellchat.params = list(
#         nboot = 1000,
#         exclude_anns = NULL,
#         thresh = 1,
#         assay = "RNA",
#         .normalize = FALSE,
#         .do_parallel = FALSE,
#         .raw_use = TRUE
#     ))


saveRDS(cmbc_cc, "data/output/cbmc_run.rds")
