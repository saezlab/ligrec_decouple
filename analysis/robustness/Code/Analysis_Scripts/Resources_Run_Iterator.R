#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------

# This is the implementation of Run_Iterator.R that runs with slurm inputs.



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
  require(magrittr)


  # Define the functions needed to perform our analysis

  # Define the iterator wrapper function (wrap_resource_Iterator), which
  # produces the end results. The wrapper iterates the evaluator function and
  # collates its results.
  source("analysis/robustness/Code/Resource_Dilution/RD_Iterator.R")
  # Define general functions for data processing in the iterator
  source("analysis/robustness/Code/Resource_Dilution/Iterator_Processing_Functions.R")
  # Define functions for capturing metadata and saving iterator results
  source("analysis/robustness/Code/Resource_Dilution/Iterator_Metadata_and_Saves.R")


  # The evaluator function tests method robustness over multiple dilution stages
  # The script below defines resource_Robustness, a function for testing
  # resource robustness
  source("analysis/robustness/Code/Resource_Dilution/RD_Evaluator.R")
  # Define functions that help the evaluator run LIANA
  source("analysis/robustness/Code/Resource_Dilution/Evaluator_LIANA_Functions.R")


  # Define the dilutor function, a function that dilutes resources.
  source("analysis/robustness/Code/Resource_Dilution/RD_Dilutor.R")
  # Define general helper functions for the dilutor
  source("analysis/robustness/Code/Resource_Dilution/Dilutor_Processing_Functions.R")


  # Define utility functions that are used in resource dilution and in cluster
  # reshuffling.
  # Define common functions for the Iterators
  source("analysis/robustness/Code/Utilities/Iterator_Functions.R")
  # Define common functions for console outputs and plots
  source("analysis/robustness/Code/Utilities/User_Outputs_and_Plots.R")
  # Define common functions that help with NATMI parameters
  source("analysis/robustness/Code/Utilities/LIANA_Utilities.R")

}



#------------------------------------------------------------------------------#
# 2. Get Resource Dilution Robustness Results ----------------------------------

# We get SLURM input arguments, including job_id to mark NATMI files with.
args = commandArgs(TRUE)

modify_baseline <- as.logical(args[1])
job_id <- as.numeric(args[2])


# First we load testdata from the data folder.
# We also give a label (testdata_type, choose "seurat_pbmc" or "liana_test")
testdata_type <- "seurat_pbmc"
testdata <- extract_Testdata(testdata_type = testdata_type)


# Double-check the Inputs
{
  important_args <-
    c("testdata_type",
      "modify_baseline",
      "job_id")

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
robustness_resource <-
  wrap_resource_Iterator(
    testdata = testdata,
    testdata_type = testdata_type,
    modify_baseline = modify_baseline,
    NATMI_tag = job_id
  )


# Extra info for console output
print(robustness_resource)
