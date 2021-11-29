# ----------------------------------------------------------------------------#
# A. Overview ------------------------------------------------------------------

# This script serves the purpose of retrieving all the data you need to perform 
# the analyses in this repository. It retrieves, unpacks and then stores the 
# data in the Data folder. At the moment, only one external data source is used.

# ----------------------------------------------------------------------------#
# B. Loading Packages ----------------------------------------------------------

require(tidyverse)


# ----------------------------------------------------------------------------#
# C. PBMC Data -----------------------------------------------------------------

# Setting up the PBMC data set download (url of data and location to download)
# Make sure your working directory is the ligrec_robustness folder

url_PBMC <- 
  "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"

path_PBMC <- 
  str_glue(getwd(), "/Data/pbmc3k_filtered_gene_bc_matrices.tar.gz")

# PBMC data is downloaded to the data folder of ligrec_robustness
download.file(url_PBMC, path_PBMC)

# The tarball is extracted in the data folder
untar(path_PBMC, exdir = str_glue(getwd(), "/Data"))
