# Load libs
require(tidyverse)
require(magrittr)
require(Seurat)
require(liana)
require(SPOTlight)

# brain analysis directory
brain_dir <- "data/input/spatial/brain_cortex"


# I) Prep Allen Brain Atlas
# think whether I want to use this, possibly better for LIANA, but for some reason deconv. topics are slightly less specific
cortex_sc <- readRDS(file.path(brain_dir, "allen_cortex.rds")) # whole atlas
# SPOTlight tutorial
# cortex_sc <- readRDS(glue::glue("{system.file(package = 'SPOTlight')}/allen_cortex_dwn.rds"))

# Downsample (since SPOTlight shows stable performance with 100 cells,
# we downsample the dataset for computational speed and memory)
cortex_sc@meta.data %<>%
    mutate(subclass = str_replace_all(subclass, "[/]", ".")) %>%
    mutate(subclass = str_replace_all(subclass, " ", "."))
cortex_sc <- subset(cortex_sc, cells = rownames(cortex_sc@meta.data))
Idents(cortex_sc) <- cortex_sc@meta.data$subclass
gc()

# Normalize appropriately
cortex_sc %<>%
    Seurat::SCTransform(verbose = FALSE,
                        conserve.memory= TRUE) %>%
    Seurat::RunPCA(verbose = FALSE) %>%
    Seurat::RunUMAP(dims = 1:30, verbose = FALSE)

Seurat::DimPlot(cortex_sc,
                group.by = "subclass",
                label = TRUE) + Seurat::NoLegend()

# Remove things that we don't need
cortex_sc@assays$RNA@scale.data <- matrix()
# save the formatted object
saveRDS(cortex_sc, file.path(brain_dir, "allen_cortex_prep.rds"))

# save a downsampled object (to be used for deconv?)
meta_dwn <- cortex_sc@meta.data %>%
  rownames_to_column(var = "barcode") %>%
  group_by(subclass) %>%
  slice_sample(n=200) %>%
  ungroup() %>%
  as.data.frame() %>%
  column_to_rownames("barcode")
cortex_sc <- subset(cortex_sc, cells = rownames(meta_dwn))
Idents(cortex_sc) <- meta_dwn$subclass
saveRDS(cortex_sc, file.path(brain_dir, "allen_cortex_dwn.rds"))


# save markers
cluster_markers_all <- Seurat::FindAllMarkers(object = cortex_sc,
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = cluster_markers_all,
        file = file.path(brain_dir, "markers.rds"))


# II) Deconvolute Mouse slides ----
# Read SC reference and cluster markers
cortex_sc <- readRDS(file.path(brain_dir, "allen_cortex_dwn.rds"))
cluster_markers_all <- readRDS(file.path(brain_dir, "markers.rds"))


if (! "stxBrain" %in% SeuratData::AvailableData()[, "Dataset"]) {
    # If dataset not downloaded proceed to download it
    SeuratData::InstallData("stxBrain")
}

# Deconvolute slides /w SPOTlight
slides <- c("anterior1",
            "anterior2",
            "posterior1",
            "posterior2")

slides %>%
    map(function(slide){
        message(slide)

        seurat_object <- SeuratData::LoadData("stxBrain", type = slide)
        print(seurat_object)

        set.seed(1234)
        # Run deconvolution as from tutorial and paper
        spotlight_ls <- spotlight_deconvolution(
            se_sc = cortex_sc,
            counts_spatial = seurat_object@assays$Spatial@counts,
            clust_vr = "subclass", # cell-type annotations
            cluster_markers = cluster_markers_all, # df with marker genes
            cl_n = 200, # number of cells per cell type
            hvg = 3000, # Number of HVG
            ntop = NULL, # How many of the marker genes to use (by default all)
            transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorization and NLS
            method = "nsNMF", # Factorization method
            min_cont = 0, # Remove those cells contributing to a spot below a certain threshold
            assay = "SCT"
            )
        gc()
        saveRDS(spotlight_ls, str_glue("{brain_dir}/{slide}_doconvolution2.RDS"))
    })

# check deconv nmf topics/results
deconv_results <- slides %>%
  map(function(slide){
    deconv_res <- readRDS(str_glue("{brain_dir}/{slide}_doconvolution2.RDS"))
    return(deconv_res)
  }) %>%
  setNames(slides)
nmf_mod <- deconv_results$posterior1[[1]]

h <- NMF::coef(nmf_mod[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(
  h = h,
  train_cell_clust = nmf_mod[[2]])

# Check topics/cell type specificity
topic_profile_plts[[2]]




