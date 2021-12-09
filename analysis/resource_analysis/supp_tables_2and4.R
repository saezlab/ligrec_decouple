# Load required packages
require(tidyverse)
require(OmnipathR)
require(magrittr)
require(proxy)
require(viridis)
require(RCurl)
require(liana)
require(shadowtext)
require(logger)
require(rlang)
library(patchwork)
library(ggplotify)
require(UpSetR)

# Source required functions
source("src/plot_utils.R")
source("analysis/resource_analysis/resource_descriptive.R")

# Obtain list with CCC Resources
ligrec <- readRDS("data/input/ligrec.RDS")

# Supp. Table 2 ----
# Obtain Curated interactions from resources which are fully manually curated
omni_curated <- OmnipathR::curated_ligand_receptor_interactions(
    curated_resources = c("Guide2Pharma", "ICELLNET", # "HPMR", "Kirouac2010" # Too old
                          "CellTalkDB", "CellChatDB", "connectomeDB2020"
                          ),
    cellphonedb = TRUE,
    cellinker = TRUE,
    talklr = FALSE,
    signalink = TRUE)

# Decomplexify and keep only relevant columns
omni_curated %<>%
    liana:::decomplexify() %>%
    select(source, target) %>%
    mutate(curated_flag = 1)

# Obtain curated percentages
curation_info <- select_resource("all") %>%
    compact() %>%
    map(function(res){
        res %>% liana:::decomplexify() %>%
            select(source, target) %>%
            left_join(omni_curated, by = c("source", "target")) %>%
            mutate(curated_flag = as.factor(replace_na(curated_flag, 0))) %>%
            mutate(total = n()) %>%
            group_by(curated_flag) %>%
            mutate(num =  n()) %>%
            ungroup() %>%
            select(curated_flag, total, num) %>%
            distinct() %>%
            rowwise() %>%
            mutate(perc = num/total * 100)
    }) %>%
    enframe(name="Resource") %>%
    unnest(value) %>%
    filter(curated_flag == 1) %>%
    select(Resource, Percentage = perc) %>%
    arrange(desc(Percentage))

curation_info %>%
    write.csv("tables/supp_table2_perc.csv", row.names = FALSE)

# Supp. Table 4 ----

# Format ligrec
ligrec_olap <- ligrec %>%
    ligrec_overheats %>%
    ligrec_overlap

# Get localisations
ligrec_loc <- ligrec_olap %>%
    localization_ligrec_classes

# UniProt dictionary to convert to symbols
up <- UniProt.ws::UniProt.ws(taxId=9606)


# Obtain 5 random examples for each signalling type
set.seed(1)
examples <- ligrec_loc$interactions %>%
    select(
        signalling_type=location,
        transmitter=source,
        receiver=target,
        resource) %>%
    group_by(signalling_type) %>%
    sample_n(5) %>%
    ungroup()

# Get Keys
keys <- unlist(union(str_split(examples$transmitter, pattern = "_"),
                     str_split(examples$receiver, pattern = "_")))

# Get UniProt Mapping
up_dict <- UniProt.ws::select(up, keytype = c("UNIPROTKB"),
                              columns = c("GENES","REVIEWED"),
                              keys = keys) %>%
    filter(REVIEWED == "reviewed") %>%
    dplyr::select(uniprot = UNIPROTKB,
                  genesymbol = GENES) %>%
    # Keep only first gene symbol (i.e. official one)
    tibble() %>%
    mutate(genesymbol = gsub(" .*", "", genesymbol))


examples %>%
    left_join(up_dict, by=c("transmitter"="uniprot")) %>%
    rename(transmitter_symbol = genesymbol) %>%
    left_join(up_dict, by=c("receiver"="uniprot")) %>%
    rename(receiver_symbol = genesymbol) %>%
    select(
        signalling_type,
        transmitter,
        transmitter_symbol,
        receiver,
        receiver_symbol,
        resource,
    ) %T>%
    print() %>%
    # mutate(across(ends_with("type"), ~recode(.x,
    #                                         "T"="Transmembrane",
    #                                         "P"="Peripheral",
    #                                         "S"="Secreted"))) %>%
    write.csv("tables/supp_table4.csv", row.names = FALSE)



