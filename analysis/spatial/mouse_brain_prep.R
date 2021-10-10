require(tidyverse)
require(magrittr)
require(liana)
require(SPOTlight)
require(Seurat)

# I) Prep Allen Brain Atlas
cortex_sc <- readRDS("data/input/spatial/brain_cortex/allen_cortex.rds")

# Downsample (since SPOTlight shows stable performance with 100 cells,
# we downsample the dataset for computational speed and memory)
cortex_sc@meta.data %<>%
    mutate(subclass = str_replace_all(subclass, "[/]", ".")) %>%
    mutate(subclass = str_replace_all(subclass, " ", ".")) %>%
    rownames_to_column(var = "barcode") %>%
    group_by(subclass) %>%
    slice_sample(n=200) %>%
    ungroup() %>%
    as.data.frame() %>%
    column_to_rownames("barcode")
cortex_sc <- subset(cortex_sc, cells = rownames(cortex_sc@meta.data))
Idents(cortex_sc) <- cortex_sc@meta.data$subclass
gc()

# Normalize appropriately
cortex_sc %<>%
    Seurat::SCTransform(verbose = FALSE) %>%
    Seurat::RunPCA(verbose = FALSE) %>%
    Seurat::RunUMAP(dims = 1:30, verbose = FALSE)

Seurat::DimPlot(cortex_sc,
                group.by = "subclass",
                label = TRUE) + Seurat::NoLegend()

# save the downsampled and formatted object
saveRDS(cortex_sc, "data/input/spatial/brain_cortex/allen_cortex_dwn.rds")


# save markers
cluster_markers_all <- Seurat::FindAllMarkers(object = cortex_sc,
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = cluster_markers_all,
        file = "data/input/spatial/brain_cortex/markers.rds")


# II) Deconvolute Mouse slides ----
# Read SC reference and cluster markers
cortex_sc <- readRDS("data/input/spatial/brain_cortex/allen_cortex_dwn.rds")
cluster_markers_all <- readRDS("data/input/spatial/brain_cortex/markers.rds")


if (! "stxBrain" %in% SeuratData::AvailableData()[, "Dataset"]) {
    # If dataset not downloaded proceed to download it
    SeuratData::InstallData("stxBrain")
}

# Deconvolute slides /w SPOTlight
c("anterior1", "anterior2", "posterior1", "posterior2") %>%
    map(function(slide){
        message(slide)

        seurat_object <- SeuratData::LoadData("stxBrain", type = slide)
        print(seurat_object)

        # Run deconvolution as from tutorial and paper
        spotlight_ls <- spotlight_deconvolution(
            se_sc = cortex_sc,
            counts_spatial = seurat_object@assays$Spatial@counts,
            clust_vr = "subclass", # Variable in sc_seu containing the cell-type annotation
            cluster_markers = cluster_markers_all, # Dataframe with the marker genes
            cl_n = 200, # number of cells per cell type to use
            hvg = 3000, # Number of HVG to use
            ntop = NULL, # How many of the marker genes to use (by default all)
            transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
            method = "nsNMF", # Factorization method
            min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
        )

        saveRDS(spotlight_ls, str_glue("data/input/spatial/brain_cortex/{slide}_doconvolution.RDS"))
    })



saveRDS(object = spotlight_ls, file = here::here("data/spotlight_ls.rds"))




anterior <- SeuratData::LoadData("stxBrain", type = "anterior1")
anterior2 <- SeuratData::LoadData("stxBrain", type = "anterior2")
posterior <- SeuratData::LoadData("stxBrain", type = "posterior1")
posterior2 <- SeuratData::LoadData("stxBrain", type = "posterior2")



# III) Run LIANA ----
cortex_sc <- readRDS("data/input/spatial/brain_cortex/allen_cortex_dwn.rds")
murine_resource <- readRDS("data/input/murine_omnipath.RDS")

liana_res <- liana_wrap(cortex_sc,
                        cellchat.params=list(organism="mouse"),
                        resource = "custom",
                        external_resource = murine_resource,
                        assay = "SCT",
                        expr_prop=0.1)

# squidpy sets gene names to upper (in the processing), revert this to title (i.e. murine)
# set all others to title just in case
liana_res %<>% map(function(res) res %>%
                       mutate_at(.vars = c("ligand", "receptor"), str_to_title))
saveRDS(liana_res, "data/spotlight_liana.rds")
