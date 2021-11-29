#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # This is the CR_LIANA_Functions.R script. It defines functions that 
  # run LIANA in the iterator.
}



#------------------------------------------------------------------------------#
# 1. Define Functions ----------------------------------------------------------

# liana_with_warnings()
{
  #' Run liana with warning handling support
  #' 
  #' @param testdata A seurat object to run liana_wrap on.
  #' 
  #' @param methods_vector The methods liana_wrap should apply.
  #' 
  #' @param liana_warnings Either TRUE, FALSE or "divert". Should the warnings 
  #' from LIANA++, which are often repetitive and unhelpful, be either 
  #' suppressed or alternatively diverted to the Logs folder? When these types 
  #' of warning are left in, they can often displace valuable warnings. Be 
  #' careful with this setting, as suppressing warnings is obviously risky.
  #' 
  #' @param warning_logile Where should the warnings be logged? Only necessary
  #' when liana_warnings == "divert".
  #'  
  #' @param tag An additional tag that specifies the NATMI save files to this
  #' iterator run.
  #' 
  #' @param ... Variable arguments to be passed to liana_wrap.
  #' 
  #' @return A list named after methods_vector that contains a tibble of CCIs
  #' for each method.
  
  
  liana_with_warnings <- function(testdata,
                                  methods_vector,
                                  
                                  tag,
                                  
                                  liana_warnings,
                                  warning_logfile,
                                  ...) {
    
    
    # If the user wants warnings, simply run liana_wrap
    if (liana_warnings == TRUE) {
      
      # NATMI needs custom file paths, check LIANA_Utilities.R for detail.
      liana_results <-
        liana_wrap(testdata,
                   method   = methods_vector,
                   resource = c('OmniPath'),
                   call_natmi.params = unique_natmi_filepaths(tag),
                   ...)
      
      
    # If the user wants to divert the warnings to a log file, use 
    # divert_Warnings() and liana_wrap to do it.
    } else if (liana_warnings == "divert") {
      
      divert_Warnings({
        
        # NATMI needs custom file paths, check LIANA_Utilities.R for detail.
        liana_results <-
          liana_wrap(testdata,
                     method   = methods_vector,
                     resource = c('OmniPath'),
                     call_natmi.params = unique_natmi_filepaths(tag),
                     ...)
        
      }, logFile = warning_logfile)
      
      
    # If the user doesn't want warnings, run liana_wrap inside suppressWarnings
    } else if (liana_warnings == FALSE) {
      
      suppressWarnings({
        
        # NATMI needs custom file paths, check LIANA_Utilities.R for detail.
        liana_results <-
          liana_wrap(testdata,
                     method   = methods_vector,
                     resource = c('OmniPath'),
                     call_natmi.params = unique_natmi_filepaths(tag),
                     ...)
        
      })
      
      
    }
    
    
    # If only one method is fed to liana_wrap(), the output data structure is
    # different. So here we convert it to the same data structure that would
    # exist if multiple methods had been called, so the code below works
    # properly for single method runs too.
    if (length(methods_vector) == 1) {
      
      liana_results        <- list(liana_results)
      names(liana_results) <- methods_vector
      
      
    }
    
    return(liana_results)
    
  }
  
  
}


# iterate_liana_wrap()
{
  #' Run liana_with_warnings() for multiple mismatch proportions and seeds
  #' 
  #' @description This function modifies the input testdata so that it has 
  #' the cluster annotations and metadata derived from reshuffled_clusters. 
  #' It then runs LIANA++ on the modified metadata. This process is performed
  #' for every mismatch proportion and seed provided (as long as the 
  #' corresponding metadata is in reshuffled_clusters).
  #' 
  #' In the process, the testdata is subset so it only contains cells that are 
  #' labelled within its metadata. If the metadata was reshuffled, this should
  #' change nothing, but it also enables this function to work with subset 
  #' metadata.
  #' 
  #' @param testdata A seurat object to run liana_with_warnings on.
  #' 
  #' @param methods_vector The methods liana_with_warnings should apply.
  #' 
  #' @param reshuffled_clusters A list of meta.data tables from testdata that
  #' have been reshuffled to various degrees (shuffle_Clusters() outputs). 
  #' Should be categorized by mismatch proportion and seed.
  #'   
  #' @param mismatch_props As a named list of proportions between 0 and 1. To 
  #' what degree were the cluster annotations mismatched? Will be used to 
  #' iterate over and to subset reshuffled_clusters.
  #' 
  #' @param seed_list As a named list of seeds. Will be used to iterate 
  #' over and to subset reshuffled clusters.
  #' 
  #' @param liana_warnings Either TRUE, FALSE or "divert". Should the warnings 
  #' from LIANA++, which are often repetitive and unhelpful, be either 
  #' suppressed or alternatively diverted to the Logs folder? When these types 
  #' of warning are left in, they can often displace valuable warnings. Be 
  #' careful with this setting, as suppressing warnings is obviously risky.
  #' 
  #' @param warning_logile Where should the warnings be logged? Only necessary
  #' when liana_warnings == "divert".
  #' 
  #' @param tag An additional tag that specifies the NATMI save files to this
  #' iterator run.
  #' 
  #' @param ... Variable arguments to be passed to liana_wrap.
  #' 
  #' @return All the liana_results
  
  
  iterate_liana_wrap <- function(testdata,
                                 methods_vector,
                                 
                                 mismatch_props,
                                 seed_list,
                                 reshuffled_clusters,
                                 
                                 tag,
                                 
                                 liana_warnings,
                                 warning_logfile,
                                 ...) {
    
    # Start up runtime measurements.
    runtime <- list("Start Iterations" = Sys.time())
    
    # Let the user keep up in the output
    print_Title(
      "2. LIANA for Mismatched Cluster Annotations",
      space_after = 0,
      super = TRUE
    )
    
    # for every mismatch_prop and seed, run liana_with_warnings()
    liana_results <- map(mismatch_props, function(mismatch_prop) {
      
      liana_results_mismatch <- map(seed_list, function(seed) {
        
        # let the user know what parameters the liana++ outputs are for
        print_Title(str_glue(mismatch_prop * 100,
                             " % Mismatch  --  Iteration ",
                             seed,
                             ": ",
                             as.character(Sys.time())))
        
        # Print to the log file what parameters the liana++ warnings refer to
        if (liana_warnings == "divert") {
          cat("\n\n\n",
              str_glue("|=======================================",
                       "=======================================|"),
              "\n\n",
              str_glue("  ",
                       mismatch_prop * 100,
                       " % Mismatch  --  Iteration ",
                       seed,
                       ": ",
                       as.character(Sys.time())),
              "\n\n\n",
              file = warning_logfile,
              append = TRUE)
          
          
        }
        
        
        
        
        # Create anew testdata with the reshuffled meta.data
        reshuffled_testdata <- testdata
        
        mismatch_name       <- names(which(mismatch_props == mismatch_prop))
        
        reshuffled_testdata@meta.data <-
          reshuffled_clusters[[mismatch_name]][[seed]]
        
        # Make sure only cells labelled in the metadata are in the testdata
        reshuffled_testdata <- 
          subset(reshuffled_testdata, 
                 cells = rownames(reshuffled_testdata@meta.data))
        
        # Let the suer know the impact of subsetting
        cat(str_wrap(
          str_glue(
            "The testdata was subset to only contain cells that ",
            "are labelled in the modified metadata. ",
            "As a result of this, ",
            ncol(testdata) - ncol(reshuffled_testdata),
            " cells were removed."),
          80),
          "\n\n")
        
        # Set the idents of the reshuffled testdata through the new metadata
        Idents(reshuffled_testdata) <-
          reshuffled_testdata@meta.data$cluster_key
        




    
        
        # Run liana with the modified testdata
        liana_results_mismatch_seed <-
          liana_with_warnings(
            liana_warnings  = liana_warnings,
            methods_vector  = methods_vector,
            testdata        = reshuffled_testdata,
            warning_logfile = warning_logfile,
            tag = tag,
            ...
          )
        
        return(liana_results_mismatch_seed)
        
      })
      
      return(liana_results_mismatch)
      
    })
    
    # Runtime checkpoint
    runtime[["Shuffled Clusters"]] <- Sys.time()
    
    # Let the use know these outputs are going to be for default LIANA
    print_Title("3. LIANA with Default Annotations.",
                super = TRUE)
    
    # Run liana with the original testdata
    original_results <-
      liana_with_warnings(
        testdata        = testdata,
        methods_vector  = methods_vector,
        liana_warnings  = liana_warnings,
        warning_logfile = warning_logfile,
        tag = tag,
        ...
      )
    
    # Create a similar format for the original results as the reshuffled ones
    original_results <-
      map(seed_list, function(seed) {
        return(original_results)
      })
    
    # runtime checkpoints
    runtime[["Default Clusters"]] <- Sys.time()
    
    
    # format outputs
    complete_liana_results <- original_results %>%
      list("Reshuffle_0" = .)  %>%
      append(., liana_results) %>%
      append(., list("runtime" = runtime))
    
    return(complete_liana_results)
    
  }
  
  
}



