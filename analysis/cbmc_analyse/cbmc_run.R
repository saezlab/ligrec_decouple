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
    seurat_object = seurat_object,
    resource = c("OmniPath", "CellPhoneDB"),
    natmi.params = list(
        expr_file = "em.csv",
        meta_file = "metadata.csv",
        output_dir = "NATMI_cbmc",
        assay = "RNA",
        num_cor = 4,
        .format = TRUE,
        .write_data = FALSE,
        .seed = 1004,
        .natmi_path = "~/Repos/LIANA/NATMI"),
    squidpy.params = list(
        cluster_key="rna_annotations",
        n_perms=10000,
        threshold=0.1,
        seed=as.integer(1004)
    ),
    cellchat.params = list(
        nboot = 1000,
        exclude_anns = NULL,
        thresh = 1,
        assay = "RNA",
        .normalize = FALSE,
        .do_parallel = FALSE,
        .raw_use = TRUE
    ))


saveRDS(cmbc_cc, "data/output/cbmc_run.rds")
