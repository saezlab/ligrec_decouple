# Run on local

# Load and and format Data
SeuratData::InstallData("cbmc")

cbmcdata <- SeuratData::LoadData("cbmc") %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures()

cbmcdata@meta.data %<>%
    mutate(celltype = as.factor(rna_annotations)) %>%
    # remove + as it breaks Squidpy
    mutate(rna_annotations = str_replace_all(rna_annotations, "[+]", "")) %>%
    mutate(rna_annotations = str_replace(rna_annotations, "/", " ")) %>%
    mutate(rna_annotations = as.factor(rna_annotations))
cbmcdata <- subset(cbmcdata, cells = rownames(cbmcdata@meta.data))
cbmcdata <- SetIdent(cbmcdata, value = cbmcdata@meta.data$rna_annotations)
saveRDS(cbmcdata, "data/input/cbmc_seurat.rds")
cbmcdata <- readRDS("data/input/cbmc_seurat.rds")


# Call liana
cmbc_cc <- liana_wrap(seurat_object = cbmcdata,
           resource = "all",
           cellchat.params = list(
                   nboot = 1000,
                   thresh = 1,
                   assay = "RNA"
                   ),
           squidpy.params = list(
               cluster_key="seurat_annotations",
               n_perms=10000,
               threshold=0.01,
               seed=as.integer(1004)
           ))



default_args <- list(
    "cellchat" = cellchat.params %<>%
        `%||%`(list(
            nboot = 100,
            exclude_anns = NULL,
            thresh = 1,
            assay = assay,
            .normalize = FALSE,
            .do_parallel = FALSE,
            .raw_use = TRUE
        )),

    'connectome' = connectome.params %<>%
        `%||%`(list(
            min.cells.per.ident = 1,
            p.values = TRUE,
            calculate.DOR = FALSE,
            assay = assay,
            .format = TRUE,
            .spatial = FALSE
        )),

    'italk' = italk.params %<>%
        `%||%`(list(
            assay = assay,
            .format = TRUE,
            .DE = TRUE
        )),

    'natmi' = natmi.params %<>%
        `%||%`(list(
            expr_file = "em.csv",
            meta_file = "metadata.csv",
            output_dir = "NATMI_test",
            assay = "RNA",
            num_cor = 4,
            .format = TRUE,
            .write_data = TRUE,
            .seed = 1004,
            .natmi_path = NULL
        )),

    'sca' = sca.params %<>%
        `%||%`(list(
            assay = assay,
            .format = TRUE,
            s.score = 0,
            logFC = log2(1.5))),

    'squidpy' = squidpy.params %<>%
        `%||%`(list(
            cluster_key="seurat_annotations",
            n_perms=1000,
            threshold=0.01,
            seed=as.integer(1004)
        ))
)



# 1. Squidpy -------------------------------------------------------------------
squidpy_results <- call_squidpyR(seurat_object = cbmcdata,
                                 python_path = "/home/dbdimitrov/anaconda3/envs/theisverse/bin/python",
                                 omni_resources = omni_resources,
                                 .ident = "rna_annotations")
saveRDS(squidpy_results, "output/cbmc_res/squidpy_results.rds")


# 2. NATMI --------------------------------------------------------------------
# save OmniPath Resource to NATMI format
natmi_results <- call_natmi(omni_resources = omni_resources,
                            seurat_object = cbmcdata,
                            wd_path = "/home/dbdimitrov/Repos/ligrec_decoupleR",
                            omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
                            natmi_path = "~/Repos/NATMI",
                            em_path = "~/Repos/ligrec_decoupleR/input/cbmcdata_em.csv",
                            ann_path = "~/Repos/ligrec_decoupleR/input/cbmcdata_ann.csv",
                            output_path = "~/Repos/ligrec_decoupleR/output/cbmc_res/natmi_results",
                            .write_data = TRUE,
                            .subsampling_pipe = FALSE,
                            .assay = "RNA",
                            .num_cor = 11
)
saveRDS(natmi_results, "output/cbmc_res/natmi_results.rds")


# 3. SCA ----------------------------------------------------------------------
sca_results <- omni_resources %>%
    map(function(db)
        call_sca(op_resource = db,
                 cbmcdata,
                 assay = 'RNA',
                 .format = TRUE,
                 s.score = 0,
                 logFC = log2(1.5)
        ))
saveRDS(sca_results, "output/cbmc_res/sca_results.rds")

# 4. Connectome ----------------------------------------------------------------
conn_results <- omni_resources %>%
    map(function(db)
        call_connectome(seurat_object = cbmcdata,
                        .spatial = FALSE,
                        op_resource = db,
                        min.cells.per.ident = 1,
                        p.values = TRUE,
                        calculate.DOR = FALSE,
                        .format = TRUE,
                        assay = 'RNA'))
saveRDS(conn_results, "output/cbmc_res/conn_results.rds")

# 5. iTALK ---------------------------------------------------------------------
italk_results <- omni_resources %>%
    map(function(db)
        call_italk(op_resource = db,
                   cbmcdata,
                   assay = 'RNA',
                   .format = TRUE,
                   .DE = TRUE
        ))
saveRDS(italk_results, "output/cbmc_res/italk_results.rds")


# 6. CellChat -----------------------------------------------------------------
cellchat_results <- omni_resources %>%
    map(function(db) call_cellchat(op_resource = db,
                                   seurat_object = cbmcdata,
                                   nboot = 1000,
                                   exclude_anns = c(),
                                   thresh = 1,
                                   assay = "RNA",
                                   .normalize = FALSE,
                                   .do_parallel = FALSE)) %>%
    setNames(names(omni_resources))
saveRDS(cellchat_results, "output/cbmc_res/cellchat_results.rds")
