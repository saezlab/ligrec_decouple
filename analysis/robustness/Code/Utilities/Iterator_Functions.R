#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # The Iterator_Functions.R script. As a utilities script, it defines functions
  # that are used in both resource dilution and cluster reshuffling. The 
  # functions here are grouped because they all play a role in running the 
  # iterator for each of these respective processes.
  
}



#------------------------------------------------------------------------------#
# 1. Define Functions ----------------------------------------------------------

# extract_Testdata()
{
  #' Helper function that gets a specific seurat object from the outputs folder
  #' 
  #' @param testdata_type As a string. Which testdata should be retrieved? 
  #' Either "seurat_pbmc" or "liana_test". Seurat_pbmc is the data set used in 
  #' the seurat tutorial, while liana_test is the testdata that comes with 
  #' LIANA++, and is a small subset of seurat_pbmc. 
  #' 
  #' @return A seurat object loaded from the outputs folder or liana package.
  
  
  extract_Testdata <- function(testdata_type) {
    
    # Get seurat or liana test data
    if (testdata_type == "seurat_pbmc") {
      
      # Read testdata from outputs
      testdata <- readRDS(file = "Data/pbmc3k_final.rds")     
      
    } else if (testdata_type == "liana_test") {
      
      # Where is the liana testdata located?
      liana_path <- system.file(package = 'liana')       
      # Read the testdata from its location.
      testdata <- 
        readRDS(file.path(liana_path, "testdata", "input", "testdata.rds"))   
      
      # removing superfluous values
      rm(liana_path)
      
      
    } else {
      
      # error if its not one of the supported data sets
      stop("Testdata name not recognized!")
      
    }
    
    
    # Return the seurat.
    return(testdata)
    
  } # end of function
}


# rank_overlap()
{
  #' Takes get_n_top_ranks outputs that have an LR_ID and determines their 
  #' overlap
  #' 
  #' @description An LR_ID is a unique identifier for a CC Interaction. It is
  #' simply the concatenated source name, target name, ligand name and receptor
  #' name. A CCI is fully characterized by its LR_ID, which is why they can be
  #' used when comparing two tibbles with top ranked interactions. 
  #'
  #' @param main_ranks A tibble of of top ranked interactions
  #' @param comparison_ranks A tibble of top ranked interactions
  #' @param verbose Should the function describe the overlap to you or not?
  #' @param expect_same_size A boolean that indicates if it is expected that the
  #' two tibbles being compared are of the same size. If this is not expected,
  #' (FALSE) it is assumed the user knows that overlaps are directional and no 
  #' warning is produced when the two tibbles are not the same size.
  #'
  #' @return The overlap (0-1) between the two input frames in contents of the
  #' LR_ID column, as well as an optional print statement that gives more 
  #' detail.
  
  rank_overlap <-
    function(main_ranks, 
             comparison_ranks, 
             verbose = TRUE, 
             expect_same_size = TRUE) {
      
      # calculate overlap between LR_IDs
      overlap <- sum(comparison_ranks$LR_ID %in% main_ranks$LR_ID)
      percentage_overlap <- overlap / nrow(main_ranks)
      
      
      # describe the output to the user
      if (verbose == TRUE) {
        print(str_glue(""))
        cat(
          str_wrap(
            str_glue(
              "The main ranking and the comparison ranking have ",
              as.character(overlap),
              " LR_IDs in common. Which corresponds to a ",
              as.character(round(percentage_overlap * 100, 2)),
              "% overlap."), 
            width = 60), "\n")
      }
      
      
      # put out a warning if the rankings are not of the same length, in this 
      # case the overlap percentage is only based on the main_ranks, which might
      # catch the user off-guard
      if ((expect_same_size == TRUE) &&
          (nrow(main_ranks) != nrow(comparison_ranks))) {
        warning("Rankings are not of same length.
            Percentage based on main_ranks argument.")
      }
      
      
      return(percentage_overlap)
      
    } #end of function
  
  
}


# calculate_Runtime()
{
  #' Converts a list of checkpoint names and associated system times into a 
  #' convenient tibble that highlights the time that passed between the
  #' checkpoints
  #' 
  #' @description Takes the runtime output of resource_Robustness(), which is a 
  #' named list of checkpoints in time, and creates a  output tibble that has 
  #' the time between each checkpoint in it and the time elapsed up until that
  #' checkpoint.
  #' 
  #' @param runtime A list of checkpoints from the Iterator and the times at 
  #' which they were reached. 
  #' 
  #' @return A tibble giving an overview of the runtime of a piece of code, 
  #' typically of the Iterator..
  
  
  calculate_Runtime <- function(runtime) {
    
    # save the names of the time-points for later
    runtime_labels <- names(runtime)
    
    # convert run time to numeric so we can perform arithmetic operations on
    # them. In this case we need it for subtractions, to calculate the duration
    # between checkpoints
    runtime_numeric <- as.numeric(runtime)
    
    # We calculate the passage of time between checkpoints in the 
    # resource_Robustness().
    # Step duration is the duration of a step between neighboring checkpoints.
    # Time elapsed is the duration between the completion of a step and the 
    # start of the script.
    
    step_duration <- c(0) # No time has passed when the script is initialized.
    time_elapsed  <- c(0) # No time has passed when the script is initialized.
    
    # starting with the second index of runtime_numeric until the last index
    for (i in 2:length(runtime_numeric)) {
      
      # subtract the preceding checkpoint from the checkpoint at i, this is the 
      # amount of time that passed between these two checkpoints
      step_duration <- c(step_duration, 
                         runtime_numeric[[i]] - runtime_numeric[[i-1]])
      
      # subtract the very first checkpoint from the checkpoint at i, this is all
      # the time that has elapsed up until now.
      time_elapsed  <- c(time_elapsed,
                         runtime_numeric[[i]] - runtime_numeric[[1]])
      
    }
    
    # Turn seconds into time periods using lubridate and round for simplicity
    # Time periods are HH:MM:SS, which is earier to understand than just values
    # in seconds.
    step_duration <- round(seconds_to_period(step_duration))
    time_elapsed  <- round(seconds_to_period(time_elapsed))
    
    
    # summarize all the runtime data in a tibble
    runtime <- runtime               %>%
      as_tibble_col()                %>%
      unnest(cols = c(value))        %>%
      rename("Start Time" = "value") %>% 
      add_column("Step Name"      = runtime_labels, .before = 1) %>%
      add_column("Step Duration"  = step_duration) %>%
      add_column("Time Elapsed"   = time_elapsed) 
    
    
    # Get rid of clutter in the environment
    rm(runtime_numeric, 
       step_duration, 
       time_elapsed, 
       runtime_labels,
       i)
    
    # return the runtime tibble
    return(runtime)
  }
}

