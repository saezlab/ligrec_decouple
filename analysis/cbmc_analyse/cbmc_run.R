# Run on local
require(liana)
require(tidyverse)
require(magrittr)

# Load and and format Data
SeuratData::InstallData("cbmc")

cbmcdata <- SeuratData::LoadData("cbmc")

cbmcdata@meta.data %<>%
    mutate(celltype = as.factor(rna_annotations)) %>%
    # remove + as it breaks Squidpy
    mutate(rna_annotations = str_replace_all(rna_annotations, "[+]", "")) %>%
    mutate(rna_annotations = str_replace(rna_annotations, "/", " "))  %>%
    mutate(rna_annotations = if_else(rna_annotations=="B", "B cell", rna_annotations)) %>%
    filter(rna_annotations != "Mouse")
cbmcdata <- subset(cbmcdata, cells = rownames(cbmcdata@meta.data))
cbmcdata <- Seurat::SetIdent(cbmcdata, value = cbmcdata@meta.data$rna_annotations)
seurat_object <- cbmcdata %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures()
seurat_object@meta.data %<>%
    mutate(rna_annotations = as.factor(rna_annotations))
saveRDS(seurat_object, "data/input/cbmc_seurat.rds")

# Call liana
cmbc_cc <- liana_wrap(
    seurat_object = cbmcdata,
    resource = c("OmniPath", "CellPhoneDB"),
    method = c("cellchat", "italk", "natmi", "sca", "squidpy"),
    natmi.params = list(
        expr_file = "em.csv",
        meta_file = "metadata.csv",
        output_dir = "NATMI_test",
        assay = "RNA",
        num_cor = 4,
        .format = TRUE,
        .write_data = FALSE,
        .seed = 1004,
        .natmi_path = "~/Repos/LIANA/NATMI"
        )
    )



# CellChat
res1 <- call_cellchat(
    op_resource = NULL,
    seurat_object = seurat_object,
    nboot = 10,
    exclude_anns = NULL,
    thresh = 1,
    assay = "RNA",
    .normalize = FALSE,
    .do_parallel = FALSE,
    .raw_use = TRUE
)

# SCA
res1 <- call_sca(op_resource = NULL,
                 seurat_object = seurat_object,
                 assay = 'RNA',
                 .format = TRUE,
                 s.score = 0,
                 logFC = log2(1.5))

# iTALK
res1 <- call_italk(op_resource = NULL,
                   seurat_object = seurat_object,
                   assay = 'RNA',
                   .format = TRUE,
                   .DE = TRUE)

# NATMI
res1 <- call_natmi(op_resource = select_resource("OmniPath"),
                   seurat_object = seurat_object,
                   expr_file = "em.csv",
                   meta_file = "metadata.csv",
                   output_dir = "NATMI_cbmc",
                   assay = "RNA",
                   num_cor = 4,
                   .format = TRUE,
                   .write_data = TRUE,
                   .seed = 1004,
                   .natmi_path = "~/Repos/LIANA/NATMI")

# squidpy
res1 <- call_squidpy(seurat_object = seurat_object,
                     op_resource = select_resource("OmniPath"),
                     cluster_key="rna_annotations",
                     n_perms=10000,
                     threshold=0.01,
                     seed=as.integer(1004))

# Connectome
res1 <- call_connectome(
    seurat_object = seurat_object,
    op_resource = select_resource("OmniPath")[[1]], # Default = No sig hits
    .spatial = FALSE,
    min.cells.per.ident = 1,
    p.values = TRUE,
    calculate.DOR = FALSE,
    assay = 'RNA',
    .format = TRUE
)


# Call liana
cmbc_cc <- liana_wrap(
    seurat_object = cbmcdata,
    resource = c("OmniPath", "CellPhoneDB"),
    method = c("connectome")
)
