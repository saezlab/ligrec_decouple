# res_liana_with_warnings()
{
  #' Run liana with warning handling support
  #' 
  #' @description This is a similar version of the cluster reshuffling function.
  #' 
  #' @param testdata A seurat object to run liana_wrap on.
  #' 
  #' @param methods_vector The methods liana_wrap should apply.
  #' 
  #' @param resource The resource parameter to be passed to liana_wrap.
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
  
  
  res_liana_with_warnings <- function(testdata,
                                      methods_vector,
                                      resource,
                                      external_resource = NULL,
                                      
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
                   resource = resource,
                   call_natmi.params = unique_natmi_filepaths(tag),
                   external_resource = external_resource,
                   ...)
      
      
      # If the user wants to divert the warnings to a log file, use 
      # divert_Warnings() and liana_wrap to do it.
    } else if (liana_warnings == "divert") {
      
      divert_Warnings({
        
        # NATMI needs custom file paths, check LIANA_Utilities.R for detail.
        liana_results <-
          liana_wrap(testdata,
                     method   = methods_vector,
                     resource = resource,
                     call_natmi.params = unique_natmi_filepaths(tag),
                     external_resource = external_resource,
                     ...)
        
      }, logFile = warning_logfile)
      
      
      # If the user doesn't want warnings, run liana_wrap inside suppressWarnings
    } else if (liana_warnings == FALSE) {
      
      suppressWarnings({
        
        # NATMI needs custom file paths, check LIANA_Utilities.R for detail.
        liana_results <-
          liana_wrap(testdata,
                     method   = methods_vector,
                     resource = resource,
                     call_natmi.params = unique_natmi_filepaths(tag),
                     external_resource = external_resource,
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