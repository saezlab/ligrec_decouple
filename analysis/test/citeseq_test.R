require(SingleCellExperiment)
require(Seurat)
require(liana)
require(tidyverse)
require(magrittr)
source("analysis/test/citeseq_src.R")




# make this thing into a function
# filter out any interactions that don't include the ADT-RNA intersect to save time. - tried gives different results :D
# perform it with 3 datasets + 3 intervals (0, 0.05, 0.1, 0.2)
# plot together
# analyse the proportions

# filter out control ADTs,
# look for alias genesymbols


# We Run Kendal Correlation, which supposedly uses tau-b scores (i.e. accounts for ties)
# Only 7 receptors match to OmniPath LRs, if I consider also ligands it would be a higher number,
# but very similar, but it would need to be separately, i.e. by source and target, so it becomes confusing
# and in my opinion redundant - ligand and receptors overlap also.




# Load matrix and convert to Seurat object ----
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.1.0/5k_pbmc_protein_v3
mats <- Seurat::Read10X("data/input/10x_citeseq/5k_pbmc_protein_v3_filtered_feature_bc_matrix/filtered_feature_bc_matrix/")
rownames(mats$`Antibody Capture`) %<>%
    gsub("\\_.*","", .)

seurat_object <- SeuratObject::CreateSeuratObject(counts = mats$`Gene Expression`)
seurat_object[["ADT"]] <- CreateAssayObject(counts = mats$`Antibody Capture`)

# run def analysis
seurat_object %<>%
    FindVariableFeatures() %>%
    NormalizeData() %>%
    ScaleData() %>%
    RunPCA(verbose = TRUE) %>%
    FindNeighbors(reduction = "pca") %>%
    FindClusters(resolution = 0.4, verbose = TRUE)

clust.anns <- str_glue("c{levels(seurat_object)}")
names(clust.anns) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, clust.anns)
seurat_object@meta.data <- seurat_object@meta.data %>%
    mutate(seurat_clusters = as.factor(str_glue("c{seurat_clusters}")))
Idents(seurat_object)

# Normalize ADT
seurat_object <- NormalizeData(seurat_object,
                               assay = "ADT",
                               normalization.method = "CLR")
# seurat_object <- ScaleData(seurat_object, assay = "ADT")


# save test seurat_objcet
saveRDS(seurat_object, "data/input/cmbc_seurat_test.RDS")



# RUN LIANA on cbmc test data ----
liana_res <- liana_wrap(seurat_object,
                        squidpy.params=list(cluster_key = "seurat_clusters",
                                            seed = as.integer(1)),
                        expr_prop = 0.2)
liana_res %<>% liana_aggregate()
saveRDS(liana_res, "data/output/test_citeseq_02prop.RDS")


# Run CiteSeq Correlation pipeline ----
## Load Seurat Object
seurat_object <- readRDS("data/input/citeseq/5k_pbmcs_nextgem/5k_pbmcs_nextgem_seurat.RDS")
clust.anns <- str_glue("c{levels(seurat_object)}")
names(clust.anns) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, clust.anns)
seurat_object@meta.data <- seurat_object@meta.data %>%
    mutate(seurat_clusters = as.factor(str_glue("c{seurat_clusters}")))
Idents(seurat_object)
liana_wrap <- liana_wrap(seurat_object,
                         method="call_natmi")



# Get Symbols of Receptors from OP
op_resource <- select_resource("OmniPath")[[1]]
receptor_syms <- c(op_resource$target_genesymbol)
cluster_key <- "seurat_clusters"
liana_res <- readRDS("data/input/citeseq/5k_pbmcs_nextgem/5k_pbmcs_nextgem-liana_res-0.RDS")

adt_lrcorr_summed <- wrap_adt_corr(seurat_object = seurat_object,
                                   liana_res = liana_res,
                                   op_resource = op_resource,
                                   cluster_key = "seurat_clusters")

adt_lrcorr_summed %>%
    # remove Mean/Median ranks
    filter(!(method %in% c("median_rank", "mean_rank"))) %>%
    # rename methods
    mutate(method = str_to_title(method)) %>%
    mutate(method = gsub("\\..*","", method)) %>%
    # rename metrics
    mutate(metric = if_else(metric=="scale", "Cluster-specific Mean", metric)) %>%
    mutate(metric = if_else(metric=="mean", "Mean", metric)) %>%
    mutate(metric = if_else(metric=="prop", "Cell Proportion", metric)) %>%
    ggplot(aes(x = metric, y = estimate, colour = method, shape = significant)) +
    geom_point(size = 6, position = position_dodge(w = 0.05)) +
    xlab("Expression Type") +
    ylab("Kendal's tau Correlation Coefficients") +
    theme_minimal(base_size = 24) +
    labs(colour = "Method", shape="Significant")





### Step-by-Step ----
# convert to singlecell object
sce <- SingleCellExperiment::SingleCellExperiment(
    assays=list(
        counts = GetAssayData(seurat_object, assay = "ADT", slot = "counts"),
        data = GetAssayData(seurat_object, assay = "ADT", slot = "data")
        ),
    colData=DataFrame(label=seurat_object@meta.data[[cluster_key]])
)

# Filter ADT Controls
sce <- sce[!(rownames(sce) %in% c("IgG1", "IgG2a", "IgG2b")),]



# Get gene aliases ----
# Obtain adt genesymbols
adt_aliases <- get_adt_aliases(adt_symbols = rownames(sce))
adt_aliases


# Get ADT stats and bind to LIANA res ----
adt_means <- get_adt_means(sce = sce,
                           receptor_syms = op_resource$target_genesymbol,
                           adt_aliases)



# Correlation between ADT-Gene (for receptors as baseline) ----
assays_corr <- get_assays_corr(seurat_object,
                               cluster_key = "seurat_clusters",
                               adt_means = adt_means)




# Get LR-ADT Correlations ----
# bind ADT means to LIANA res
liana_adt <- liana_format_adt(liana_res = liana_res,
                              adt_means = adt_means)

set.seed(2) # actually does nothing
adt_corr <- get_adt_correlations(liana_adt)


#
adt_corr %<>% bind_rows(assays_corr) %>%
    mutate(metric = gsub("\\_.*","", metric))

# plot
adt_corr %>%
    # remove Mean/Median ranks
    filter(!(method %in% c("median_rank", "mean_rank"))) %>%
    # rename methods
    mutate(method = str_to_title(method)) %>%
    mutate(method = gsub("\\..*","", method)) %>%
    # rename metrics
    mutate(metric = if_else(metric=="scale", "Cluster-specific Mean", metric)) %>%
    mutate(metric = if_else(metric=="mean", "Mean", metric)) %>%
    mutate(metric = if_else(metric=="prop", "Cell Proportion", metric)) %>%
    ggplot(aes(x = metric, y = estimate, colour = method, shape = significant)) +
    geom_point(size = 6, position = position_dodge(w = 0.05)) +
    xlab("Expression Type") +
    ylab("Kendal's tau Correlation Coefficients") +
    theme_minimal(base_size = 24) +
    labs(colour = "Method", shape="Significant")



# Get RNA-ADT corr ----
# get adt aliases
sce_rna <- SingleCellExperiment::SingleCellExperiment(
    assays=list(
        counts = GetAssayData(seurat_object, assay = "RNA", slot = "counts"),
        data = GetAssayData(seurat_object, assay = "RNA", slot = "data")
    ),
    colData=DataFrame(label=seurat_object@meta.data[[cluster_key]])
)

# called adt_means but simply runs means on sce object and formats the output
mean_summary <- scuttle::summarizeAssayByGroup(sce_rna,
                                               ids = colLabels(sce_rna),
                                               assay.type = "data",
                                               statistics = c("mean"))
rna_means <- mean_summary@assays@data$mean %>%
    # gene mean across cell types
    as_tibble(rownames = "entity_symbol") %>%
    mutate(entity_symbol = stringr::str_to_upper(entity_symbol)) %>%
    # Join alias names and rename entity symbol to newly-obtained alias symbol
    # purpose being to be able to join to the gene symbols in the RNA assay
    # left_join(adt_aliases, by = c("entity_symbol"="adt_symbol")) %>%
    # dplyr::select(-c(entity_symbol)) %>%
    # keep only receptors present in OP
    filter(entity_symbol %in% receptor_syms) %>%
    pivot_longer(-entity_symbol,
                 names_to = "target",
                 values_to = "rna_mean") %>%
    group_by(entity_symbol) %>%
    # obtain across cluster scaled means
    mutate(rna_scale = scale(rna_mean)) %>%
    unnest(rna_scale) %>%
    ungroup()



rna_adt_join <- adt_aliases %>%
    left_join(adt_means, by = c("alias_symbol"="entity_symbol")) %>%
    left_join(rna_means, by = c("alias_symbol" = "entity_symbol", "target")) %>%
    na.omit()


mean_corr <- cor.test(rna_adt_join$adt_mean,
                      rna_adt_join$rna_mean,
                      method = "kendal") %>%
    broom::tidy() %>%
    mutate(metric = "mean")

scale_corr <- cor.test(rna_adt_join$adt_scale,
                       rna_adt_join$rna_scale,
                       method = "kendal") %>%
    broom::tidy() %>%
    mutate(metric = "scale")






bind_rows(mean_corr,
          scale_corr)
# Get Summary per clust and Join alias_symbols ----
test_summ <- scuttle::summarizeAssayByGroup(sce,
                                            ids = colLabels(sce),
                                            assay.type = "data",
                                            statistics = c("mean",
                                                           "prop.detected"))



test_summ@colData
means <- test_summ@assays@data$mean %>% # gene mean across cell types
    as_tibble(rownames = "entity_symbol") %>%
    mutate(entity_symbol = stringr::str_to_upper(entity_symbol)) %>%
    # Join alias names and rename entity symbol to newly-obtained alias symbol
    # purpose being to be able to join to the gene symbols in the RNA assay
    left_join(adt_aliases, by = c("entity_symbol"="adt_symbol")) %>%
    dplyr::select(-c(entity_symbol)) %>%
    dplyr::select(entity_symbol = alias_symbol, everything())
means

props <- test_summ@assays@data$prop.detected %>% # gene prop
    as_tibble(rownames = "entity_symbol") %>%
    mutate(entity_symbol = stringr::str_to_upper(entity_symbol)) %>%
    left_join(adt_aliases, by = c("entity_symbol"="adt_symbol")) %>%
    dplyr::select(-c(entity_symbol)) %>%
    dplyr::select(entity_symbol = alias_symbol, everything())
props




# Format Means and Proportion Stats
# props
adt_props_entity <- props %>%
    filter(entity_symbol %in% receptor_syms) %>%
    pivot_longer(-entity_symbol,
                 names_to = "cluster",
                 values_to = "adt_prop")

# means
adt_means_entity <- means %>%
    filter(entity_symbol %in% receptor_syms) %>%
    pivot_longer(-entity_symbol,
                 names_to = "cluster",
                 values_to = "adt_mean") %>%
    group_by(entity_symbol) %>%
    mutate(adt_scale = scale(adt_mean)) %>%
    unnest(adt_scale) %>%
    ungroup()

# join props to means
adt_means_entity %<>%
    left_join(adt_props_entity) %>%
    dplyr::rename("target" = cluster)


# Basic correlations between LRs and ADT means -----
liana_res2 <- liana_res %>%
    filter(receptor %in% adt_means_entity$entity_symbol) %>%
    dplyr::select(source, ligand, target, receptor, ends_with("rank")) %>%
    pivot_longer(c(ligand,receptor),
                 names_to = "type",
                 values_to = "entity_symbol") %>%
    filter(type == "receptor") %>%
    distinct()

liana_adt <- liana_res2 %>%
    left_join(adt_means_entity) %>%
    dplyr::select(ends_with("rank"), starts_with("adt")) %>%
    pivot_longer(-c(adt_scale, adt_mean, adt_prop))


# run corrs -----
set.seed(1)
adt_corr <- liana_adt %>%
    group_by(name) %>%
    nest() %>%
    mutate(mean_model = map(data, ~corr_fun(.x, var="adt_mean"))) %>%
    mutate(scale_model = map(data, ~corr_fun(.x, var="adt_scale"))) %>%
    mutate(prop_model = map(data, ~corr_fun(.x, var="adt_prop")))


# Check correlations
mean_corr <- adt_corr %>%
    dplyr::select(name, mean_model) %>%
    unnest(cols = c(mean_model)) %>%
    mutate(mean_mlog10p = -log10(p.value)) %>%
    dplyr::select(name, mean_pval = p.value, mean_estimate = estimate)

scale_corr <- adt_corr %>%
    dplyr::select(name, scale_model) %>%
    unnest(cols = c(scale_model)) %>%
    dplyr::select(name, scale_pval = p.value, scale_estimate = estimate)

prop_corr <- adt_corr %>%
    dplyr::select(name, prop_model) %>%
    unnest(cols = c(prop_model)) %>%
    mutate(prop_mlog10p = -log10(p.value)) %>%
    dplyr::select(name, prop_pval = p.value, prop_estimate = estimate)

corr_format <- mean_corr %>%
    # join
    left_join(scale_corr) %>%
    left_join(prop_corr) %>%
    dplyr::rename("method" = name) %>%
    pivot_longer(-method,
                 names_to = "corr_type") %>%
    separate(corr_type, into = c("metric", "stat")) %>%
    pivot_wider(names_from = stat,
                values_from = value) %>%
    mutate(significant = if_else(pval <= 0.05, TRUE, FALSE)) %>%
    dplyr::select(-pval) %>%
    ungroup()

# plot
corr_format %>%
    # remove Mean/Median ranks
    filter(!(method %in% c("median_rank", "mean_rank"))) %>%
    # rename methods
    mutate(method = str_to_title(method)) %>%
    mutate(method = gsub("\\..*","",method)) %>%
    mutate(estimate = -1*estimate) %>% # reverse negative estimates
    # rename metrics
    mutate(metric = if_else(metric=="scale", "Cluster-specific Mean", metric)) %>%
    mutate(metric = if_else(metric=="mean", "Mean", metric)) %>%
    mutate(metric = if_else(metric=="prop", "Cell Proportion", metric)) %>%
    ggplot(aes(x = metric, y = estimate, colour = method, shape = significant)) +
    geom_point(size = 6, position = position_dodge(w = 0.05)) +
    xlab("Expression Type") +
    ylab("Kendal's tau Correlation Coefficients") +
    theme_minimal(base_size = 24) +
    labs(colour = "Method", shape="Significant")


# ROC according to Proportion thresholds ----
# hist(adt_means_entity[,-1] %>% as.matrix() %>% as.numeric(), breaks=100)
# normalize and cut anything below 0, i.e. assume that its noise

adt_props_entity
liana_prop <- liana_res2 %>%
    left_join(adt_props_entity, by = c("target" = "cluster", "entity_symbol"))



### Run on Murine -----
seurat_object <- readRDS("data/input/citeseq/spleen_lymph_206//spleen_lymph_206_seurat.RDS")
op_resource <- liana::select_resource("OmniPath")[[1]] %>%
    convert_to_murine()

murine_liana_206 <- liana_wrap(seurat_object,
                               resource = "custom",
                               external_resource = op_resource,
                               expr_prop = 0.1,
                               cellchat.params=list(organism="mouse"))

# everything to title (some that are nottotitle will be mismatched otherwise, due to the upper of squidpy)
murine_liana %<>% map(function(res) res %>% mutate_at(.vars = c("ligand", "receptor"), str_to_title))

murine_liana %<>% liana_aggregate
saveRDS(murine_liana_206, "data/input/citeseq/spleen_lymph_206/spleen_lymph_206-liana_res-0.1.RDS")



# 111 ----
seurat_object <- readRDS("data/input/citeseq/spleen_lymph_101/spleen_lymph_101_seurat.RDS")
op_resource <- liana::select_resource("OmniPath")[[1]] %>%
    convert_to_murine()

murine_liana_111 <- liana_wrap(seurat_object,
                               resource = "custom",
                               external_resource = op_resource,
                               expr_prop = 0.1,
                               cellchat.params=list(organism="mouse"))

murine_liana_111$squidpy %<>%
    mutate_at(.vars = c("ligand", "receptor"), str_to_title)
murine_liana_111 %<>% liana_aggregate
saveRDS(murine_liana_111, "data/input/citeseq/spleen_lymph_101//spleen_lymph_111-liana_res-0.1.RDS")


xd <- load_adt_lr(dir = citeseq_dir,
                  subdir = "spleen_lymph_101",
                  op_resource = murine_resource,
                  cluster_key = "seurat_clusters",
                  liana_pattern = "liana_res-0.1.RDS",
                  organism = "mouse"
                  )



### Manually convert
genes <- rownames(seurat_object)[rownames(seurat_object) %in%
                            union(op_resource$target_genesymbol,
                                  op_resource$source_genesymbol)]



prots <- rownames(seurat_object@assays$ADT) %>%
    as_tibble() %>%
    filter(!str_detect(value, "Ctrl"))  %>%
    filter(!str_detect(value, "Ligand")) %>%
    mutate(adt = gsub("^[^-]+-\\s*", "", value)) %>%
    mutate(adt = gsub("-[^-]*$", "", adt)) %>%
    mutate(adt = gsub(")", "", adt)) %>%
    separate(adt, into = c("adt", "add"), sep = "\\(")

syms <- get_adt_aliases(str_to_title(prots$adt),
                        organism = "human") %>%
    bind_rows(get_adt_aliases(str_to_title(prots$adt),
                              organism = "mouse")) %>%
    mutate(across(everything(), str_to_title)) %>%
    distinct()


murine_adt <- get_alias_table(organism = "mouse")

syms[str_to_title(syms$alias_symbol) %in% str_to_title(genes),]
syms[str_to_title(syms$adt_symbol) %in% str_to_title(genes),]


genes[order(genes)]

# Manual annotations obtained via NCBI, ENSMBL, and from the ADT's descriptions
# With priority given to the ADT description
# Cd107a = "Lamp1"
# Cd11a = "Itgal"
# Cd11b = "Itgam"
# Cd11c = "Itgax"
# Cd120b = "Tnfrsf1b"
# Cd137 = "Tnfrsf9"
# Cd14 = "Cd14"
# Cd140a = "Pdgfra"
# Cd152 = Ctla4
# Cd159a = "Klrc1"
# Cd16-32 = "Fcgr3"
# Cd160 = "By55"
# Cd163 = "Cd163"
# Cd169 = c("Siglec-1", "Siglec1", "Sn")
# Cd170 = c("Singlec-F")
# Cd172 = "Sirpa"
# Cd182 = "Cxcr2"
# Cd183 = "Cxcr3"
# Cd185 = "Cxcr5"
# Cd186 = "Cxcr6"
# Cd198 = "Cxcr8"
# Cd200R3 = "Cd200r3"
# Cd201 = "Procr"
# Cd204 = "Msr1"
# Cd207 = "CD207"
# Cd21-Cd35 = c("Cr2","Cr1")
# Cd22 = "Cd22"
# Cd226 = "Dnam-1"
# Cd253 = "Trail"
# Cd26 = c("DPP-4","Dpp4", "Dpp-4)
# Cd270 = c("Hvem", "Tnfrsf14")#
# Cd272 = c("Btla")
# Cd273 = "Pdcd1lg2"
# Cd274 = "Cd274"
# Cd278 = Cd278
# Cd279 = c("B7h1", "Cd279", "Pdcd1l1", "Pdcd1lg1", "Pdl1")
# Cd28 = "Cd28"
# Cd300lg = Cd300lg
# Cd301a = c("Mgl", "Mgl1", "CD301", "CD301a")
# Cd301b = c("Mgl2", "CD301b")
# Cd304 = c("Nrp1")
# Cd309 = c("CD309", "Vegfr2-Flk-1", "Kdr")
# Cd31 = c("Pecam1", "Pecam-1")

prots <- rownames(readRDS("data/input/citeseq/adt_raw.RDS")) %>%
    as_tibble() %>%
    filter(!str_detect(value, "Ctrl"))  %>%
    filter(!str_detect(value, "Ligand")) %>%
    mutate(adt = gsub("^[^_]+_\\s*", "", value)) %>%
    mutate(adt = gsub("_[^_]*$", "", adt)) %>%
    mutate(adt = gsub(")", "", adt)) %>%
    separate(adt, into = c("adt", "add"), sep = "\\(")

manual_list <- list(
    CD105 = c("Eng", "Endo"),
    CD107a = c("Lamp1", "Perk"),
    CD120b = c("Tnfrsf1b", "Tnfr2"),
    CD16.32 = "Fcgr3",
    CD198 = "Cxcr8",
    CD199 = c("Ccr9", "Cmkbr10"),
    CD201 = c("Procr", "Epcr"),
    CD21.CD35 = c("Cr2","Cr1"),
    CD278.1 = c("Icos, Ailim"),
    CD300c.d = c("Cd300c", "Clm6"),
    CD301a = c("Mgl", "Mgl1"),
    CD301b= c("Mgl2"),
    CD309.1 = c("Kdr","Flk1"),
    CD326 = c("Epcam", "Tacstd1"),
    CD34.1 = "CD34",
    CD45.1 = c("Ptprc"),
    CD45.2 = c("Ptprc"),
    CD45R.B220 = c("Ptprc"),
    CD49a = c("Itga1"),
    CD90.1 = c("Thy1", "Thy-1"),
    CD90.2 = c("Thy1", "Thy-1"),
    D62E = c("Sele", "Elam-1"),
    F4.80 = c("Adgre4", "Emr4"),
    FceRIa = c("Fcer1g", "Fce1g"),
    FolateReceptorb = c("Folr2", "Fbp2", "Folbp2"),
    IL.21Receptor = c("Il21r", "Nilr"),
    IL.33Ra = c("Il1rl1", "Ly84", "St2", "Ste2"),
    Ly.49A = c("Klra1", "Ly-49", "Ly-49a", "Ly49", "Ly49A"),
    Ly.6A.E = c("Sca-1", "Ly6e"),
    Ly.6G = "Ly6g6e",
    Ly.6G.Ly.6C = "Gr-1",
    MAdCAM.1 = "Madcam1",
    Mac.2 = "Lgals3",
    NK.1.1 = c("Klrb1c", "Ly55c", "Nkrp1c"),
    P2X7R = "P2rx7",
    # leave out T-cell Receptor variants
    # TCRVb13.1 = "Trbv13",
    # TCRVb5.1,
    # TCRVb8.1,
    # TCRVr1.1.Cr4,
    # TCRVr2,
    # TCRVr3,
    # TCRbchain,
    # TCRr.d,
    Tim.4 = c("Timd4", "Tim4"),
    anti.P2RY12 = c("P2ry12", "P2y12")#,
    # integrinb7,
    # B6,
    # B6.1
)

adt_mat <- GetAssayData(seurat_object, assay = "ADT", slot = "data")



adt_mat <- adt_mat[!str_detect(rownames(adt_mat), "Ctrl"),]

