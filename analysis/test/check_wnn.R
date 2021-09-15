require(Seurat)
require(SingleCellExperiment)

require(tidyverse)
require(liana)

op_resource <- select_resource("OmniPath")[[1]]

# WNN
## Check metadata
wnn_meta <- readRDS("data/input/seurat_wnn_data/metadata.rds") # Seurat4 is required
unique(wnn_meta$celltype.l1)
unique(wnn_meta$celltype.l2) # 31 - too many. preferably exclude some?
# e.g. rm eryth, Doublet, all naive cell types (i.e. - 3);  merge cDC 1 and 2, Doublet,
unique(wnn_meta$celltype.l3) # way too many

# 228 antibodies (ADTs)
wnn_meta <- readRDS("data/input/seurat_wnn_data/metadata.rds")
nrow(wnn_meta)
assay <- readRDS("data/input/seurat_wnn_data/adt_expression.rds")



sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts = assay),
                                                  colData=DataFrame(label=wnn_meta$celltype.l2))

test_summ <- scuttle::summarizeAssayByGroup(sce,
                                            ids = colLabels(sce),
                                            assay.type = "counts")




test_summ@colData
means <- test_summ@assays@data$mean %>% # gene mean across cell types
    as_tibble(rownames = "entity_symbol")
means

props <- test_summ@assays@data$prop.detected %>% # gene prop
    as_tibble(rownames = "entity_symbol")
props

op_syms <- c(op_resource$source_genesymbol,
             op_resource$target_genesymbol)

adt_props_entity <- props %>%
    filter(entity_symbol %in% op_syms)

adt_means_entity <- means %>%
    filter(entity_symbol %in% op_syms)

# looks okay I guess, obviously need to be normalized before that
adt_scaled_entity <- adt_means_entity %>%
    pivot_longer(-entity_symbol) %>%
    group_by(entity_symbol) %>%
    mutate(value = scale(value)) %>%
    unnest(value) %>%
    pivot_wider(names_from = name,
                values_from = value)


# adt_ligs <- means %>%
#     filter(entity_symbol %in% op_resource$source_genesymbol)
#
# adt_recepts <- means %>%
#     filter(entity_symbol %in% op_resource$target_genesymbol)
#
