#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # Welcome! This is the Run_Iterator.R script for the cluster reshuffling 
  # process. The goal of this script is to measure the robustness of CCI 
  # (Cell-Cell Interaction) inference methods when the cluster annotations of 
  # your data become more and more inaccurate. scRNA CCI inference methods 
  # assume that the cell clustering is 100 % accurate, which of course it never 
  # really can be.
  # 
  # The way this robustness measurement is achieved in the following steps:
  # 
  # 1. Take the existing cluster annotations and reshuffle them to various
  # degrees of mismatch. Since standard reshuffling can lead to a cluster
  # annotation replaced by itself, and we want to exactly quantify the 
  # annotation-mismatch to prediction-mismatch relationship, we use a special 
  # form of reshuffling. When we reshuffle x % this special way, the mismatch 
  # is also exactly x %. Since reshuffling has inherent randomness to it, we
  # produce multiple reshuffled annotations for every mismatch percentage. 
  # This lets us see the spread of possible predictions based on that degree
  # of mismatch.
  # 
  # 2. Run all the methods in LIANA++ on testdata with the standard cluster 
  # annotations. This will be our baseline of CCI predictions. We then run
  # LIANA++ on testdata with reshuffled cluster annotations. These will be our
  # comparative values.
  # 
  # 3. Each method has its own scoring system and measures of significance. To
  # unify their rankings we simply take the highest ranked interactions for
  # each method. These are considered top-ranked interactions. 
  # 
  # 4. Once we have the top-ranked interactions for each method and degree of
  # cluster annotation mismatch, we calculate the overlap between the baseline
  # top-ranked CCI predictions and the predictions with reshuffled annotations.
  # This tells us what percentage of of the original CCI's predicted are in the 
  # reshuffled CCI predictions. 
  # 
  # It is notable that these top-ranked interactions being compared  aren't 
  # always the same size. For one, we make sure to include items at the cutoff 
  # point that are tied, and for another some methods have inbuilt significance
  # filters. These won't put out any interactions once a certain mismatch 
  # proportion in the annotations is reached, because none of them meet the 
  # method's standards for significance. 
  # 
  # Since overlap is normed by the size of one but not both of the objects 
  # involved, this could in theory pose a problem to accuracy. If for example 
  # the reshuffled predictions were somehow much larger due to ties they may 
  # recreate the original predictions, along with a bunch of faulty ones, which
  # cannot be considered accurate.  However, in our case the comparative CCIs 
  # are always of similar size or smaller than the original predictions. When 
  # it comes to their ability to recreate their original predictions, these 
  # types of size differences don't negatively impact the analysis. 
  # 
  # 5. Once the overlap to the original predictions for each reshuffling
  # proportion is known it is summarized in one tibble, plotted, and saved to
  # the outputs folder.
  # 
  # For more insight on the programmatic approach, please read the 
  # CR_Iterator.R script.
  
}



#------------------------------------------------------------------------------#
# 1. Setup ---------------------------------------------------------------------
{
  # In this segment, we set up all the required infrastructure to run the
  # wrapper. The default parameters for the wrapper are set in
  # CR_Iterator.R, but custom parameters can be set in the function call.
  
  # Load required Packages
  require(tidyverse)
  require(Seurat)
  require(liana)
  require(lubridate)
  
  
  # Define the functions needed to perform our analysis
  
  # Define wrap_cluster_Iterator(), a function that performs the core analysis
  source("Code/Cluster_Reshuffling/CR_Iterator.R")
  # Define Iterator helper functions that work with the top-ranked CCI's
  source("Code/Cluster_Reshuffling/Iterator_Top_Ranks.R")
  # Define Iterator helpers that work at a meta level, saving results and so on 
  source("Code/Cluster_Reshuffling/Iterator_Meta_and_Saves.R")

  
  # Define functions for generating the many reshuffled cluster annotations
  source("Code/Cluster_Reshuffling/CR_Shuffler_Functions.R")
  # Define functions for subsetting cluster annotations, instead of reshuffling
  source("Code/Cluster_Reshuffling/CR_Subsetter_Functions.R")
  # Define functions for iterating LIANA on all the reshuffled clusters
  source("Code/Cluster_Reshuffling/CR_LIANA_Functions.R")


  # Define utility functions that help the iterator run 
  source("Code/Utilities/Iterator_Functions.R")
  # Define utiltiy functions for user-facing console outputs and plotting
  source("Code/Utilities/User_Outputs_and_Plots.R")
  # Define common functions that help with NATMI parameters
  source("Code/Utilities/LIANA_Utilities.R")
  
  
}



#------------------------------------------------------------------------------#
# 2. Get Cluster Reshuffling Robustness Results --------------------------------

# First we load testdata from the data folder. 
# We also give a label (testdata_type, choose "seurat_pbmc" or "liana_test")
testdata_type <- "liana_test"  
testdata      <- extract_Testdata(testdata_type = testdata_type)

# We need a tag to distinctly mark this analysis from others in the NATMI save
# files, such as a SLURM JOBID. If you're only running one analysis at a time
# you don't need this tag
tag <- 12345

# We run the wrapper with test settings
robustness_reshuffle_default <- 
  wrap_cluster_Iterator(testdata       = testdata,
                        testdata_type  = testdata_type,
                        NATMI_tag      = tag,
                        methods_vector = c("call_sca", 
                                           "call_natmi",
                                           "call_connectome"),
                        number_seeds   = 2,
                        top_n          = 100,
                        mismatch_props = c(0.2, 0.4, 0.6),
                        trial_run      = TRUE)
  
  



