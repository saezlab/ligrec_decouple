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

# Obtain Curated interactions from resources which are fully manually curated
omni_curated <- OmnipathR::curated_ligand_receptor_interactions(
    curated_resources = c("Guide2Pharma", "ICELLNET", "HPMR", "Kirouac2010", # Too old
                          "CellTalkDB", "CellChatDB", "connectomeDB2020"
    ),
    cellphonedb = TRUE,
    cellinker = TRUE,
    talklr = TRUE,
    signalink = TRUE)

# Decomplexify and keep only relevant columns
omni_curated %<>%
    liana:::decomplexify() %>%
    select(source, target) %>%
    mutate(curated_flag = 1)


### Supp Figure 2 Panel D ----
# keep only the curated interactions for each resource
ligrec_curated <- ligrec %>%
    transpose() %>%
    pluck("interactions") %>%
    map(function(res)
        res %>%
            left_join(omni_curated, by=c("source", "target")) %>%
            mutate(curated_flag = as.factor(curated_flag)) %>%
            filter(curated_flag == 1) %>%
            select(source, target) %>%
            distinct()
        )

# Reproduce Supp Fig 2A., but with curated alone
binarize_curated <- ligrec_curated %>%
    binarize_resources(c("source", "target"))

# Save Figure
overheat_save(interactions_shared(binarize_curated),
              file.path("figures", str_glue("curated_shared_heat.pdf")),
              str_glue("Contained % \nCurated Interactions"))



# Supp. Table 2 (Overlap) ----
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
    arrange(Resource)

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




# Response to R3A2 ----
# Obtain Curated interactions from resources which are fully manually curated
omni_curated <- OmnipathR::curated_ligand_receptor_interactions(
    curated_resources = c("Guide2Pharma", "ICELLNET", "HPMR", "Kirouac2010", # Too old
                          "CellTalkDB", "CellChatDB", "connectomeDB2020"
    ),
    cellphonedb = TRUE,
    cellinker = TRUE,
    talklr = TRUE,
    signalink = TRUE)

# Get all interactions by DB
sep_cols <- c(str_glue("source{rep(1:50)}"))
omni_curated_db <- omni_curated %>%
    select(source, target, sources) %>%
    separate(sources, sep=";",
             into = sep_cols) %>%
    pivot_longer(-c(source, target),
                 names_to = "source_number",
                 values_to = "db") %>%
    na.omit() %>%
    select(source, target, db)


# Dependency dictionary
dict_list <- list(
    "Baccin2019" = c("Ramilowski2015", "KEGG", "Reactome"),
    "CellCall" = c("connectomeDB2020", "Cellinker", "CellTalkDB", "CellChatDB", "STRING"),
    "CellChatDB" = c("KEGG"),
    "Cellinker" = c("CellPhoneDB", "Guide2Pharma", "HMPR", "DLRP"),
    "CellPhoneDB" = c("Guide2Pharma", "I2D", "IntAct", "UniProt", "HPIDB"),
    "CellTalkDB" = c("STRING"),
    "connectomeDB2020" = c("Ramilowski2015", "CellPhoneDB", "Baccin2019", "LRdb", "ICELLNET"),
    "EMBRACE" = c("Ramilowski2015"),
    "Guide2Pharma" = c("Literature"),
    "HPMR" = c("Literature"),
    "ICELLNET" = c("STRING", "Ingenuity", "BioGRID", "Reactome", "CellPhoneDB"),
    "iTALK" = c("Ramilowski2015", "Guide2Pharma", "HMPR"),
    "Kirouac2010" = c("COPE"),
    "LRdb"= c("Ramilowski2015", "Guide2Pharma", "HPMR", "HPRD", "Reactome", "UniProt"),
    "Ramilowski2015" = c("DLRP", "HMPR", "IUPHAR", "HPRD", "STRING")#,
    # "scConnect" = c("Guide2Pharma")
    )

dict_list[[1]]

# Obtain Curated by source/db
curated_by_db <- omni_curated_db %>%
    filter(str_detect(db , dict_list[[1]])) %>%
    mutate(curated_flag = 1)

#
ligrec$Baccin2019$interactions %>%
    liana:::decomplexify() %>%
    select(source, target) %>%
    left_join(curated_by_db) %>%
    mutate(curated_flag = as.factor(replace_na(curated_flag, 0))) %>%
    mutate(total = n()) %>%
    group_by(curated_flag, db) %>%
    mutate(num =  n()) %>%
    summarise(perc = num/total * 100) %>%
    distinct() %>%
    ungroup() %>%
    mutate(db=replace_na(db, "Unknown"))

xd <- imap(dict_list, function(dbs, resource){
    curated_by_db <- omni_curated_db %>%
        filter(str_detect(db , dict_list[[resource]])) %>%
        mutate(curated_flag = 1)

    ligrec[[resource]]$interactions %>%
        liana:::decomplexify() %>%
        select(source, target) %>%
        left_join(curated_by_db, by = c("source", "target")) %>%
        mutate(curated_flag = as.factor(replace_na(curated_flag, 0))) %>%
        mutate(total = n()) %>%
        group_by(db) %>%
        mutate(num =  n()) %>%
        summarise(perc = num/total * 100) %>%
        distinct() %>%
        ungroup() %>%
        mutate(db=replace_na(db, "Own curation / Unknown")) %>%
        mutate(resource = resource)
    }) %>% bind_rows()

xd_wide <- xd %>% pivot_wider(names_from = "resource",
                   values_from = "perc",
                   id_cols = "db"
                   )

xd %>% ComplexHeatmap::Heatmap()


# Revisions2 Extend Organs + Tissue ----
ligrec <- readRDS("data/input/ligrec.RDS")
upset_args <- list()

# Get Formatted ligrec
ligrec %<>%
    ligrec_decomplexify %T>%
    ligrec_overheats %>%
    ligrec_overlap %T>%
    uniq_per_res %T>%
    ligand_receptor_upset(upset_args = upset_args)


# Top 50 Organ
organ_top50 <- ligand_receptor_classes(ligrec,
                                       resource = "HPA_tissue",
                                       organ,
                                       largest = 50,
                                       filter_annot = (level %in% c("Medium", "High") & pathology=="False" & !(status %in% c(NA, "Uncertain"))),
                                       label_annot = function(x){str_to_title(str_sub(x, start = 0, end = 30))}
)

classes_enrich(organ_top50$interactions, "interactions",
               resource="HPA_organ_top50", organ, largest = 50)



# Top 50 Tissue
tissue_top50 <- ligand_receptor_classes(ligrec,
                                        resource = "HPA_tissue",
                                        tissue,
                                        largest = 50,
                                        filter_annot = (level %in% c("Medium", "High") & pathology=="False" & !(status %in% c(NA, "Uncertain")) & (tissue!="glandular cells")),
                                        label_annot = function(x){ str_split(x, "⊎", simplify = TRUE) %>% as.data.frame() %>% str_glue_data("{V2} ({V1})") %>% str_to_title()}
                              )

classes_enrich(tissue_top50$interactions, "interactions",
               resource="HPA_tissue_top50", tissue, largest = 50)

