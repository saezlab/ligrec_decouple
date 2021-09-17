require(Seurat)


# Load matrix and convert to Seurat object
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
