#### Run LIANA with all original methods (with defaults args) + multiple resources

args = commandArgs(trailingOnly=TRUE)


# Get from commandline args
brca_subtype = args[[1]] # "TNBC", "ER", "HER2"

message(str_glue("Now Running: LIANA with {brca_subtype}"))


deconv_directory <- file.path(brca_dir,
                              "deconv",
                              str_glue("{brca_subtype}_celltype_minor"))

# Load BRCA object from Spatial Directory
seurat_object <- readRDS(file.path(deconv_directory,
                                   str_glue("{brca_subtype}_celltype_minor_seurat.RDS")))


# Here, we run LIANA but with multiple resources
liana_res <- liana_wrap(seurat_object,
                        method = c('call_natmi', 'call_connectome', 'logfc',
                                   'cellchat', 'call_sca', 'squidpy', "cytotalk"),
                        resource = c("ICELLNET", "OmniPath", "CellChatDB", "CellTalkDB"),
                        # this is passed only to squidpy, cellchat, cytotalk, and logfc
                        expr_prop=0.1,
                        cellchat.params=list(nboot=1000,
                                             expr_prop = 0)) # by default as in CellChat

# save LIANA
saveRDS(liana_res,
        file.path("data/output/comparison_out/",
                  str_glue("BRCA_{brca_subtype}_liana_res.RDS")))


