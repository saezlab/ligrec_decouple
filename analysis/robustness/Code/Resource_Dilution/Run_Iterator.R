#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # Welcome! This is the Run_Iterator.R script.
  
  # This script runs a wrapper function that analyses the robustness of 
  # predictions made by SC transcriptomic CCI inference methods in relation to
  # resource change. If you want to understand the wrapper, go to 
  # RD_Robustness_Iterator.R, where it is defined. In short, it analyses 
  # robustness in the the following steps:
  
  # 1. Run every Cell-Cell-Interaction (CCI) prediction method in LIANA++ while
  #    using OmniPath (OP) and a user specified test data.
  
  # 2. Determine the highest ranked CCI's for each method
  #    -> these are the top ranks in the generic scenario ("0% dilution")
  
  # 3. Modify ("dilute") OP by replacing a proportion of its interactions with
  #    false interactions. These following rules apply for the dilution:
  
  #   -> Interactions that were top-ranked by any method cannot be replaced
  #   -> New interactions are created from genes found in the test data, and
  #      only replace interactions relevant to the test data
  #   -> All new interactions are unique, the diluted resource has no duplicates
  #   -> The diluted interactions don't include any genes that are communication
  #      partner to themselves.
  #   -> Depending on the the selected parameters, it can also be arranged that
  #      the topology of the resource is semi-preserved.
  
  #    In general, this portion has lots of adjustable arguments. What quality
  #    of genes should dilution occur with? Should the output be logged? Etc.
  
  # 4. Once we have OP diluted to various percentages, for example 10 %, 40 %
  #    and 90 %, we rerun all the LIANA++ methods for each proportion of
  #    dilution.
  
  # 5. We determine the top-ranked interactions for these dilution stages, and
  #    then compare the overlap between the new predictions and the old 0 %
  #    prediction. We also determine the rate of diluted interactions in the
  #    top_ranks, this is usually equal to 1 - the overlap.
  
  # 6. Since dilution has elements of randomness, we repeat the above process
  #    a few times so we can collate information from multiple samples. 1-5 are
  #    all handled by resource_Robustness(), a function we iterate here and that
  #    can be found in RD_Robustness_Evaluator.R.
  
  # 7. We plot the overlap over the dilution proportion for each method,
  #    combining all the iterations into one box plot. Then we save the plots
  #    and all the results to automatically generated, descriptive file names.
  
}



#------------------------------------------------------------------------------#
# 1. Setup ---------------------------------------------------------------------
{
  # In this segment, we set up all the required infrastructure to run the 
  # wrapper. The default parameters for the wrapper are set in 
  # RD_Robustness_Iterator.R, but custom parameters can be set in the function 
  # call.
  
  # Load required Packages
  require(tidyverse)
  require(Seurat)
  require(liana)
  require(lubridate)
  
  
  # Define the functions needed to perform our analysis
  
  # Define the iterator wrapper function (wrap_resource_Iterator), which 
  # produces the end results. The wrapper iterates the evaluator function and 
  # collates its results.
  source("Code/Resource_Dilution/RD_Iterator.R")
  # Define general functions for data processing in the iterator
  source("Code/Resource_Dilution/Iterator_Processing_Functions.R")
  # Define functions for capturing metadata and saving iterator results
  source("Code/Resource_Dilution/Iterator_Metadata_and_Saves.R")
  
  
  # The evaluator function tests method robustness over multiple dilution stages
  # The script below defines resource_Robustness, a function for testing 
  # resource robustness
  source("Code/Resource_Dilution/RD_Evaluator.R")
  # Define functions that help the evaluator analyse top-ranked CCIs
  source("Code/Resource_Dilution/Evaluator_Top_Ranks_Functions.R")
  # Define functions that help the evaluator run LIANA
  source("Code/Resource_Dilution/Evaluator_LIANA_Functions.R")
  
  
  # Define the dilutor function, a function that dilutes resources.
  source("Code/Resource_Dilution/RD_Dilutor.R")
  # Define general helper functions for the dilutor 
  source("Code/Resource_Dilution/Dilutor_Processing_Functions.R")
  
  
  # Define utility functions that are used in resource dilution and in cluster
  # reshuffling.
  # Define common functions for the Iterators
  source("Code/Utilities/Iterator_Functions.R")
  # Define common functions for console outputs and plots 
  source("Code/Utilities/User_Outputs_and_Plots.R")
  # Define common functions that help with NATMI parameters
  source("Code/Utilities/LIANA_Utilities.R")
  
  
  
  
}

#------------------------------------------------------------------------------#
# 2. Get Resource Dilution Robustness Results ----------------------------------

# First we load testdata from the data folder. 
# We also give a label (testdata_type, choose "seurat_pbmc" or "liana_test")
testdata_type <- "liana_test"  
testdata      <- extract_Testdata(testdata_type = testdata_type)

# We need a tag to distinctly mark this analysis from others in the NATMI save
# files, such as a SLURM JOBID. If you're only running one analysis at a time
# you don't need this tag
tag <- 12345

# We run the wrapper with test settings
robustness_default <- 
  wrap_resource_Iterator(testdata       = testdata,
                         testdata_type  = testdata_type,
                         NATMI_tag      = tag,
                         methods_vector = c("call_connectome",
                                            "call_italk",
                                            "call_sca", 
                                            "cellchat"),
                         
                         cellchat_nperms = 10,
                         number_seeds    = 5,
                         top_n           = 500,
                         dilution_props  = c(0.2, 0.4),
                         trial_run       = TRUE,
                         modify_baseline = TRUE)

