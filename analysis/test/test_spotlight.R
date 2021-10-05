library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)

require(Seurat)
require(liana)
##

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


### LIANA on cortex -----
require(liana)
require(tidyverse)
require(magrittr)
require("biomaRt")

#' Basic function to convert human to mouse genesymbols (temporary solution)
#' @param op_resource omnipath_resource as obtained via `liana::select_resource`
#'
#' @details adapted from https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
convert_to_murine <- function(op_resource){

    # query biomaRt databases
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    # obtain tibble with human and murine genesymbol
    symbols_tibble <- getLDS(attributes = c("hgnc_symbol"),
                             filters = "hgnc_symbol",
                             values = union(op_resource$source_genesymbol,
                                            op_resource$target_genesymbol),
                             mart = human,
                             martL = mouse,
                             attributesL = c("mgi_symbol")) %>%
        dplyr::rename(human_symbol = HGNC.symbol,
                      murine_symbol = MGI.symbol) %>%
        as_tibble()

    # intentionally we introduce duplicates, if needed
    # these should be resolved when LIANA is called
    # as inappropriately matched genes will not be assigned any values
    op_resource %>%
        left_join(symbols_tibble, by=c("target_genesymbol"="human_symbol")) %>%
        mutate(target_genesymbol = murine_symbol, .keep = "unused") %>%
        left_join(symbols_tibble, by=c("source_genesymbol"="human_symbol")) %>%
        mutate(source_genesymbol = murine_symbol, .keep = "unused") %>%
        filter(!is.na(target_genesymbol) | !is.na(source_genesymbol)) %>%
        filter(!is.na(target_genesymbol)) %>%
        filter(!is.na(source_genesymbol))
}

# Run LIANA ----
# Format for LIANA
cortex_sc@meta.data %<>%
    mutate(subclass = str_replace_all(subclass, "[/]", ".")) %>%
    mutate(subclass = as.factor(subclass))
Idents(cortex_sc) <- cortex_sc@meta.data$subclass

cortex_sc <- subset(cortex_sc, cells = rownames(cortex_sc@meta.data))
Seurat::DefaultAssay(cortex_sc) <- "RNA"
cortex_sc %<>% NormalizeData()
GetAssayData(cortex_sc)


# Run
op_resource <- select_resource("OmniPath")[[1]] %>%
    convert_to_murine()

liana_res <- liana_wrap(cortex_sc,
                        cellchat.params=list(organism="mouse"),
                        resource = "custom",
                        external_resource = op_resource,
                        expr_prop=0.1)

# squidpy sets gene names to upper (in the processing), revert this to title (i.e. murine)
liana_res$squidpy %<>%
    mutate_at(.vars = c("ligand", "receptor"), str_to_title)
saveRDS(liana_res, "data/spotlight_liana.rds")


# Load Liana
liana_res <- readRDS("data/spotlight_liana.rds")
liana_res %<>% liana_aggregate()



# Spatial-Corr ROC here ----




# Simple correlation between z-transformed scores and the cell-cell proportion correlations -----
