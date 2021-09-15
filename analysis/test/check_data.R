# Single Cell ----
require(tidyverse)
require(Seurat)
require(SeuratDisk)

# CRC Korean - 18 Clusters
crc_seurat <- readRDS("data/input/crc_data/crc_korean_form.rds")
crc_seurat@meta.data

# raw CRC Kor
crc_seurat <- readRDS("data/input/crc_data/crc_korean.rds")
unique(crc_seurat@meta.data$Cell_subtype)


# + panc8 data, other SC datasets are also available singlecell package (SCE)
# NovoSparc?


# SeuratObject Cite-seq ----
cbmc_seurat <- readRDS("data/input/cbmc_seurat.rds")
cbmc_seurat@meta.data$celltype

cbmc_seurat@assays$ADT # 10 ADT features only

# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/malt_10k_protein_v3
# 18 ATDs

# https://www.10xgenomics.com/resources/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-with-cell-surface-proteins-next-gem-3-1-standard-3-1-0
# 31 surface antibodies

# IgA pilot data - comparable quality, already annotated, would Kramann agree?

# 2 datasets: 1 with ~110 antibodies and another with 200+:
# https://www.nature.com/articles/s41592-020-01050-x#data-availability



# Spatial ----
# CRC visium
crc_visium <- readRDS("data/input/crc_data/crc_visium.rds")
crc_visium@meta.data


# BC data from Rico's Collab :D
bc_visium <- readRDS("data/input/bc_data/qc_se_breast_1_full.rds")
Seurat::Idents(bc_visium)
meta <- bc_visium@meta.data %>% as_tibble()
# only proportions - no labels


# BC atlas from Paper
unique(liana_res$source)


# check disgenet
library(devtools)
install_bitbucket("ibi_group/disgenet2r") # not how to use this thing



disg <- read_delim("data/input/curated_gene_disease_associations.tsv",
           delim = "\t")

disg %>%
    filter(diseaseType == "group") %>%
    pluck("diseaseName") %>%
    unique()

# MERFISH




