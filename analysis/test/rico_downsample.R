require(dplyr)
require(tibble)
require(Seurat)
require(magrittr)
require(liana)

cluster_key <- "annotation"
rico_dir <- "/net/data.isilon/ag-saez/bq_shared/liana_debug_MI/"

# Input
seurat_object <- readRDS(file.path(rico_dir, "liana_endo_simple_obj.rds"))
# reset annotation levels
seurat_object@meta.data$annotation <- as.factor(as.character(seurat_object@meta.data$annotation))
print(seurat_object@meta.data$annotation)

# Subsample
seurat_object@meta.data %<>%
    rownames_to_column(var = "barcode") %>%
    as_tibble() %>%
    group_by(.data[[cluster_key]]) %>%
    slice_sample(prop = 0.5) %>%
    mutate( {{ cluster_key }} := as.factor(.data[[cluster_key]])) %>%
    ungroup() %>%
    as.data.frame()
rownames(seurat_object@meta.data) <- seurat_object@meta.data$barcode

# Filter
seurat_object <- subset(seurat_object, cells = rownames(seurat_object@meta.data),
                        subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 10)

# check and replace NAs
print(seurat_object@assays$RNA@data[is.na(rowSums(seurat_object@assays$RNA@data)),])
seurat_object@assays$RNA@data[is.na(seurat_object@assays$RNA@data)] <- 0


# Save object
# saveRDS(seurat_object, file.path(rico_dir, "rico_subsample1234.RDS"))

# test liana as is
liana_res <- liana_wrap(seurat_object)
# liana_res <- liana_pipe(seurat_object = seurat_object,
#                         op_resource = select_resource("OmniPath")[[1]] %>%
#                             liana:::decomplexify()
#                         )
saveRDS(liana_res, file.path(rico_dir, "rico_liana1234.RDS"))


