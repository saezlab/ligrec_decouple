require(liana)

liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

# SCA
res1 <- call_sca(op_resource = NULL,
                 seurat_object = seurat_object,
                 assay = 'RNA',
                 .format = TRUE,
                 s.score = 0,
                 logFC = log2(1.5))

# iTALK
res1 <- call_italk(op_resource = NULL,
                   seurat_object = seurat_object,
                   assay = 'RNA',
                   .format = TRUE,
                   .DE = TRUE)

# NATMI
res1 <- call_natmi(op_resource = select_resource("OmniPath"),
                   seurat_object = seurat_object,
                   expr_file = "em.csv",
                   meta_file = "metadata.csv",
                   output_dir = "NATMI_test",
                   assay = "RNA",
                   num_cor = 4,
                   .format = TRUE,
                   .write_data = FALSE,
                   .seed = 1004,
                   .natmi_path = "~/Repos/LIANA/NATMI")

# squidpy
res1 <- call_squidpy(seurat_object = seurat_object,
                     op_resource = select_resource("OmniPath"),
                     cluster_key="seurat_annotations",
                     n_perms=100,
                     threshold=0.01,
                     seed=as.integer(1004))

# Connectome
res1 <- call_connectome(
    seurat_object = seurat_object,
    op_resource = select_resource("OmniPath")[[1]], # Default = No sig hits
    .spatial = FALSE,
    min.cells.per.ident = 1,
    p.values = TRUE,
    calculate.DOR = FALSE,
    assay = 'RNA',
    .format = TRUE
)

res1 <- liana_wrap(seurat_object,
                   method = c('italk', 'sca','connectome'),
                   resource = c('OmniPath'))


