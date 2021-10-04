# check prot-rna dataset
joint <- read.csv("~/Downloads/joint.csv", header = TRUE)
unique(joint$X)

uniprots <- joint$X

liana::select_resource("OmniPath")[[1]] %>%
    select(source, source_genesymbol,
           target, target_genesymbol) %>%
    filter(source %in% uniprots) %>%
    filter(target %in% uniprots)


##
emat <- Matrix::readMM('data/input/bc_atlas/GSE176078_Wu_etal_2021_BRCA_scRNASeq/Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx')
ncol(emat)
rm(emat)
genes <- read.csv('data/input/bc_atlas/GSE176078_Wu_etal_2021_BRCA_scRNASeq/Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv')



