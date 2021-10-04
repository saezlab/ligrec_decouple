library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)

# load sc
path_to_data <- system.file(package = "SPOTlight")
cortex_sc <- readRDS(glue::glue("{path_to_data}/allen_cortex_dwn.rds"))

# load spatial
if (! "stxBrain" %in% SeuratData::AvailableData()[, "Dataset"]) {
    # If dataset not downloaded proceed to download it
    SeuratData::InstallData("stxBrain")
}

# Load data
anterior <- SeuratData::LoadData("stxBrain", type = "anterior1")

set.seed(123)
cortex_sc <- Seurat::SCTransform(cortex_sc, verbose = FALSE) %>%
    Seurat::RunPCA(., verbose = FALSE) %>%
    Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)

Seurat::DimPlot(cortex_sc,
                group.by = "subclass",
                label = TRUE) + Seurat::NoLegend()


Seurat::Idents(object = cortex_sc) <- cortex_sc@meta.data$subclass
cluster_markers_all <- Seurat::FindAllMarkers(object = cortex_sc,
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)

saveRDS(object = cluster_markers_all,
        file = "data/spotligh_test.RDS")

cluster_markers_all <- readRDS("data/spotligh_test.RDS")

set.seed(123)

spotlight_ls <- spotlight_deconvolution(
    se_sc = cortex_sc,
    counts_spatial = anterior@assays$Spatial@counts,
    clust_vr = "subclass", # Variable in sc_seu containing the cell-type annotation
    cluster_markers = cluster_markers_all, # Dataframe with the marker genes
    cl_n = 100, # number of cells per cell type to use
    hvg = 3000, # Number of HVG to use
    ntop = NULL, # How many of the marker genes to use (by default all)
    transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
    method = "nsNMF", # Factorization method
    min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
    )

saveRDS(object = spotlight_ls, file = here::here("data/spotlight_ls.rds"))
