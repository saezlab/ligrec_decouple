#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------

# This is the implementation of Run_Iterator.R that runs with slurm inputs.



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

# We get SLURM input arguments, including job_id to mark NATMI files with.
args = commandArgs(TRUE)

reshuffle_or_subset = as.character(args[1])
top_n               = as.numeric(args[2])
job_id              = as.numeric(args[3])



# First we load testdata from the data folder. 
# We also give a label (testdata_type, choose "seurat_pbmc" or "liana_test")
testdata_type <- "seurat_pbmc"  
testdata      <- extract_Testdata(testdata_type = testdata_type)


# Double-check the Inputs
{
  important_args <- 
    c("testdata_type", "reshuffle_or_subset", "top_n", "job_id") 
  
  map(important_args, function(arg) {
    c(arg, get(arg), typeof(get(arg)))
    
  }) %>% 
    transpose %>%
    set_names(c("name", "value", "data_type")) %>%
    as_tibble() %>% 
    unnest(c("name", "value", "data_type")) %>%
    print()
  
  cat("\n\nTestdata: \n")
  print(testdata)
}


# We run the wrapper with test settings
robustness_cluster <- 
  wrap_cluster_Iterator(testdata      = testdata,
                        testdata_type = testdata_type,
                        
                        reshuffle_or_subset = reshuffle_or_subset,
                        top_n               = top_n,
                        NATMI_tag           = job_id)
  

# Extra info for console output
print(robustness_cluster)

