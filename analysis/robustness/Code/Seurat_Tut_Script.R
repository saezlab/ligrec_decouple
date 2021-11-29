#------------------------------------------------------------------------------#
# A. Overview ------------------------------------------------------------------

# This a truncated version of all the code from the Seurat tutorial. The visual
# outputs have been removed, but all the processing steps are identical. If you 
# want to properly understand these processing steps and decisions, I recommend 
# reading the tutorial side by side with this code. The tutorial can be found
# here: 
#   "https://satijalab.org/seurat/articles/pbmc3k_tutorial.html"




#------------------------------------------------------------------------------#
# B. Setup ---------------------------------------------------------------------

# Loading required packages
library(tidyverse)
library(Seurat)
library(patchwork)


# Load the PBMC dataset

pbmc.data <- Read10X(data.dir = "Data/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc3k", 
                           min.cells = 3, 
                           min.features = 200)




#------------------------------------------------------------------------------#
# C. Quality Control -----------------------------------------------------------

# Adding the percentage of mitochondrial reads to metadata
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Now we subset the data using the three main QC metrics
# We cut out outliers on the high side of mt proportion and feature counts
# the plots. As recommended, we subset in one step here considering all variables.
pbmc <- 
  subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)




#------------------------------------------------------------------------------#
# D. Normalization -------------------------------------------------------------

# We normalize the counts of each cell with the total number of counts for that
# cell and multiply it by 10000.
# Additionally, we then lognormalize the entirety of the count data.
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)




#------------------------------------------------------------------------------#
# E. Variable Feature Selection ------------------------------------------------

# The most variable features are selected from the data set
# This method accounts for the fact that higher mean gene expression leads to
# higher variance by modeling the relationship of both in the data.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)




#------------------------------------------------------------------------------#
# F. Data Scaling --------------------------------------------------------------

# We adjust the gene expression such that each gene's mean expression is 0 and 
# variance is 1. This makes all genes equally weighted in analysis. Highly 
# expressed genes have as much impact as lowly expressed ones.

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)




#------------------------------------------------------------------------------#
# G. Linear Dimensional Reduction ----------------------------------------------

# PCA is performed on scaled data, in subsequent steps the number of PCs will
# be limited to a sensible estimate
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))




#------------------------------------------------------------------------------#
# H. Dimensionality Assessment -------------------------------------------------

# PCAs are constructed for multiple resamplings of the data set. The PCs that 
# most consistently show up are chosen as meta features to summarize the data 
# with. 

# NOTE: This process can take a long time for big data sets. More approximate 
# techniques such as those implemented in ElbowPlot(pbmc) can be used to reduce 
# computation time.
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)




#------------------------------------------------------------------------------#
# I. Cell Clustering -----------------------------------------------------------

# A graph-based clustering algorithm clusters the cells in the data using 10 PCs
# These cluster can be used for further analysis
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)




#------------------------------------------------------------------------------#
# J. Non-Linear Dimensional Reduction ------------------------------------------

# Non-Linear Dimensional Reduction such as UMAP can help visualize the structure
# of the data, but shouldn't be used for analysis. For completions' sake it is
# included here.

# If you haven't installed UMAP in python, you can do so via 
# reticulate::py_install(packages = 'umap-learn'). The LIANA++ conda env
# installation covers this already. 
pbmc <- RunUMAP(pbmc, dims = 1:10, 
                umap.method = 'umap-learn', 
                metric = 'correlation')


#------------------------------------------------------------------------------#
# K. Assessing Cluster Biomarkers ----------------------------------------------

# Seurat has a variety of tools for determining biomarkers, but none of these
# affect our analysis because canonical Biomarkers are known.




#------------------------------------------------------------------------------#
# L. Assigning Cluster Biomarkers ----------------------------------------------

# Canonical cluster markers are known for this data set. See the seurat tutorial
# for more information on how these labels were concluded.

new.cluster.ids <- c("Naive CD4 T", 
                     "CD14+ Mono", 
                     "Memory CD4 T", 
                     "B", 
                     "CD8 T", 
                     "FCGR3A+ Mono",
                     "NK", 
                     "DC", 
                     "Platelet")

names(new.cluster.ids) <- levels(pbmc)

pbmc <- RenameIdents(pbmc, new.cluster.ids)




#------------------------------------------------------------------------------#
# M. Saving the Seurat Object --------------------------------------------------


saveRDS(pbmc, file = str_glue(getwd(), "/Data/pbmc3k_final.rds"))


