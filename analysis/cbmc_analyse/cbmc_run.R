# Run on local
require(liana)
require(tidyverse)

# Load and and format Data
SeuratData::InstallData("cbmc")

cbmcdata <- SeuratData::LoadData("cbmc") %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures()

cbmcdata@meta.data %<>%
    mutate(celltype = as.factor(rna_annotations)) %>%
    # remove + as it breaks Squidpy
    mutate(rna_annotations = str_replace_all(rna_annotations, "[+]", "")) %>%
    mutate(rna_annotations = str_replace(rna_annotations, "/", " ")) %>%
    mutate(rna_annotations = as.factor(rna_annotations))
cbmcdata <- subset(cbmcdata, cells = rownames(cbmcdata@meta.data))
cbmcdata <- SetIdent(cbmcdata, value = cbmcdata@meta.data$rna_annotations)
saveRDS(cbmcdata, "data/input/cbmc_seurat.rds")
cbmcdata <- readRDS("data/input/cbmc_seurat.rds")

# Call liana
cmbc_cc <- liana_wrap(
    seurat_object = cbmcdata,
    resource = c("OmniPath", "CellChatDB"),
    method = "connectome",
    cellchat.params = list(
                   nboot = 10,
                   thresh = 1,
                   assay = "RNA"
                   ),
           squidpy.params = list(
               cluster_key="seurat_annotations",
               n_perms=10,
               threshold=0.01,
               seed=as.integer(1004)
               ),
           natmi.params = list(
               expr_file = "em.csv",
               meta_file = "metadata.csv",
               output_dir = "NATMI_test",
               assay = "RNA",
               num_cor = 10,
               .format = TRUE,
               .write_data = TRUE,
               .seed = 1004,
               .natmi_path = NULL
               )
           )
