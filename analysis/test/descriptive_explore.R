# Load required packages
require(tidyverse)
require(OmnipathR)
require(magrittr)
require(proxy)
require(viridis)
require(RCurl)
require(UpSetR)
require(liana)
require(shadowtext)
require(logger)
require(rlang)
library(patchwork)
library(ggplotify)

# Source required functions
source("src/plot_utils.R")
source("analysis/resource_analysis/resource_descriptive.R")

ligrec <- readRDS("data/input/ligrec.RDS")

# SignaLink check - looks good
signalink <- import_omnipath_annotations(resource = "SignaLink_pathway", wide = TRUE)
signalink_decom <- signalink %>%
    decomplexify(column = "uniprot")


# Signor
signor <- import_omnipath_annotations(resource = "SIGNOR", wide = TRUE) %>%
    decomplexify(column = "uniprot")
# ^ needs to be decomplexified!!!


# NetPath
netpath <- import_omnipath_annotations(resource = "NetPath", wide = TRUE) %>%
    decomplexify(column = "uniprot")
# ^ needs to be decomplexified!!!


# CancerSEA
cancersea <- import_omnipath_annotations(resource = "CancerSEA", wide = TRUE) %>%
    decomplexify(column = "uniprot")
# ^ needs to be decomplexified!!!


# MSigDB
msigdb <- import_omnipath_annotations(resource = "MSigDB", wide = TRUE) %>%
    decomplexify(column = "uniprot")
# ^ needs to be decomplexified!!!


# DisGeNet
disgenet <- import_omnipath_annotations(resource = "DisGeNet", wide = TRUE) %>%
    decomplexify(column = "uniprot") %>%
    # ^ needs to be decomplexified!!!
    # only curated and only somewhat disease group specific
    filter(score >= 0.3 & dpi >= 0.3 & type!="group") %>%
    filter(type!="group")


# HPA_tissue
HPA_tissue <- import_omnipath_annotations(resource = "HPA_tissue", wide = TRUE) %>%
    decomplexify(column = "uniprot")


# HGNC
HPA_tissue <- import_omnipath_annotations(resource = "HGNC", wide = TRUE) %>%
    decomplexify(column = "uniprot")


### Generate Heats
ligrec_olap <- ligrec %>%
    ligrec_decomplexify %T>%
    ligrec_overheats %>%
    ligrec_overlap

ligrec_classes_all(ligrec_olap)

