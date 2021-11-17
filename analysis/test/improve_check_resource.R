source("analysis/cytosig/cytosig_src.R")
source("src/eval_utils.R")
source("src/plot_utils.R")
source("analysis/resource_analysis/resource_descriptive.R")

require(tidyverse)
require(liana)
require(magrittr)

ligrec_olap <- ligrec %>%
    ligrec_decomplexify %T>%
    ligrec_overheats %>%
    ligrec_overlap %T>%
    uniq_per_res

op_uniques <- ligrec_olap$interactions %>%
    filter(resource == "OmniPath") %>%
    filter(unique)

op <- ligrec %>% liana:::reform_omni() %>% pluck("OmniPath")
op %>% filter(source %in% target) %>%
    filter(target %in% op_uniques$source) %>%
    print(n=100)




# Improve OmniPath ----
omni <- select_resource("OmniPath")[[1]]
flipped_interactions <- omni %>%
    filter(source %in% target) %>%
    filter((parent_intercell_source %in% c("ligand", "secreted_enzyme")))

# Get sources in flipped:
flipped_interactions$source %>% unique()
# manually checked all of these, and each was a receptor
receptors_in_ligands = c("O75462", "O95971", "P08887",
                         "P10721", "P16871", "P19438",
                         "P20333", "P22455", "P22607",
                         "P24394", "P25445", "P35916",
                         "P35968", "P42702", "Q13261")

# Get target in flipped
flipped_interactions$target %>% unique()
receptors_in_flipped <- c("P26992", "P40189", "P42702", "Q92956", "P09619",
                          "P15509", "P16144", "P19235", "P21926", "P32927",
                          "P36888", "P60033", "Q16827", "Q92729", "Q9HC73",
                          "P20333", "P19438", "P19022", "P54764", "P31785",
                          "P78552", "P00533", "P05556", "P08648", "P35968",
                          "O14786", "P05106", "P17948", "P33151", "P35916",
                          "P41231", "Q12913", "P14784", "P30530")
# ^ all of these are also receptors, since these interactions are not
# categorized as any membraine-located event -> best to filter them

# How to flip interactions:
# flip these interactions in OmniPath? ----
receptors_in_ligands

flipped_interactions %>% names

xd <- flipped_interactions %>%
    select(source, target, source_genesymbol, target_genesymbol)
xd
reflipped_interactions <- xd %>%
    # create temp names to revert
    rename_at(vars(!tidyselect::contains("target")), ~gsub("source", "target2", .x)) %>%
    rename_at(vars(!tidyselect::contains("target2")), ~gsub("target", "source2", .x)) %>%
    # revert symbol and target info
    rename_at(vars(tidyselect::contains("target2")), ~gsub("target2", "target", .x)) %>%
    rename_at(vars(tidyselect::contains("source2")), ~gsub("source2", "source", .x))

xd
reflipped_interactions %>%
    select(names(xd))

# Omni ++ (filtered) ----
omni <- select_resource("OmniPath")[[1]] %>%
    # filter mediators
    filter(!(str_detect(category_intercell_source, "cofactor")) &
               !(str_detect(category_intercell_target, "cofactor")) &
               !(str_detect(category_intercell_source, "ligand_regulator"))) %>%
    # based on the analysis above we remove ambiguous receptor-receptor interactions
    filter(!(source %in% target))

# Omni_CPDB
omni_cpdb <- select_resource("CellPhoneDB")[[1]]

# check if any ambigous interactions (wrongly annotated ligands/receptors) exist
omni_cpdb %<>%
    mutate(interaction = paste(source, target)) %>%
    mutate(interaction2 = paste(target, source)) %>%
    # OmniPath contains ~all correctly annotated receptors and ligands
    # hence we remove any ligand in CPDB that is a receptor in OmniPath
    filter(!(source %in% omni$target)) # remove mislaballed interactions (i.e duplicated ligands and receptors)

# get duplations alone
interact_dups <- c(omni_cpdb$interaction,
                   omni_cpdb$interaction2) %>%
    .[duplicated(.)]


# check if any non-membrane bound receptors/ligands are regarded
#  as both receptor and ligand
duplicated_cpdb <- omni_cpdb %>%
    filter((interaction2 %in% interact_dups) | (interaction %in% interact_dups))

# The rest, we manually fix
duplicated_cpdb$source %>% unique()
# here we include receptor proteins and such which are involved in cell-adhesion interactions
# predominantly such involved in cell-adhesion signalling (often ambiguous whether
# transmitter or receiver of signal)
actual_receivers <- c("P20292", "COMPLEX:P05556_Q13797", "P06731",
                      "P15813", "Q92478", "Q9NZS2", "P27487", "Q12918",
                      "Q9UHP7", "P41217", "Q8TD46", "Q15762", "P16150",
                       "Q15223", "Q8N126", "Q9NQS3", "P16150", "Q9BZZ2",
                      "COMPLEX:P05556_P13612"
                      )
duplicated_cpdb %>% filter(!(source %in% actual_receivers))

# Appropriately Filter CPDB ----
omni_cpdb <- select_resource("CellPhoneDB")[[1]] %>%
    # remove known receivers from source # and interactions which are duplicated but flipped
    filter(!(source %in% actual_receivers))




# remove activating_cofactor, inhibitory_cofactor

# CPDB resource ----
# Original cpdb
cpdb_complex <- read.csv("data/input/comparison/external_dbs/cpdb/complex_curated.csv")
cpdb_interaction <- read.csv("data/input/comparison/external_dbs/cpdb/interaction_curated.csv")
cpdb_hla <- read.csv("data/input/comparison/external_dbs/cpdb/hla_curated.csv")
cpdb_protein <- read.csv("data/input/comparison/external_dbs/cpdb/protein_curated.csv")

# Format stuff
cpdb_complex <- cpdb_complex %>%
    select(complex_name, uniprot_1,uniprot_2,uniprot_3, uniprot_4)  %>%
    as_tibble()

cpdb_interaction <- cpdb_interaction %>%
    select(id_cp_interaction, partner_a, partner_b, protein_name_a, protein_name_b) %>%
    as_tibble()

cpdb_protein <- cpdb_protein %>%
    select(uniprot, protein_name) %>%
    mutate(protein_name = gsub("_HUMAN", "", protein_name))

cpdb_interaction
cpdb_complex

#
orig_cpdb
