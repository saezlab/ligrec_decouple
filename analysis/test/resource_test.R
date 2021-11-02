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

# Source required functions
source("src/plot_utils.R")
source("analysis/resource_analysis/resource_descriptive.R")

## Resource descriptive analysis
# Obtain list with CCC Resources
ligrec <- compile_ligrec_descr()


# Get LIGREC olap
ligrec_olap <- ligrec %>%
    ligrec_decomplexify %T>%
    ligrec_overheats %>%
    ligrec_overlap

# Run Enrichments
ligrec_olap %T>%
    # uniq_per_res %T>%
    # ligand_receptor_upset(upset_args = list()) %>%
    ligrec_classes_all


# Check databases
OmnipathR::get_annotation_resources()

# Try DBs
TCDB <- import_omnipath_annotations(resources = "TCDB")
TopDB <- import_omnipath_annotations(resources = "TopDB")
IntOGen <- import_omnipath_annotations(resources = "IntOGen")
DGIdb <- import_omnipath_annotations(resources = "DGIdb")
HPA_subcellular <- import_omnipath_annotations(resources = "HPA_secretome")


# HPA_tissue
hpa_tissue <- import_omnipath_annotations(resources = "HPA_tissue")


hpa_tissue %>%
    filter(label=="organ") %>%
    group_by(value) %>%
    summarise(cccc = n()) %>%
    arrange(desc(cccc)) %>%
    print(n=100)


# Get DisGeNet
disgnet <- import_omnipath_annotations(resources = "DisGeNet")

xd <- disgnet %>%
    filter(label=="disease") %>%
    group_by(value) %>%
    summarise(nn = n()) %>%
    arrange(desc(nn))


msigdb <- import_omnipath_annotations(resources = "MSigDB",
                                      wide = TRUE) %>%
    filter(collection == 'hallmark')


