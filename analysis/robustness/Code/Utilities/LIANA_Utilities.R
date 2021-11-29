#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # This is the LIANA_Utilities.R script. It defines functions that 
  # help LIANA run.
}



#------------------------------------------------------------------------------#
# 1. Define Functions ----------------------------------------------------------

# unique_natmi_filepaths()
{
  #' Create unique filepaths to run with NATMI
  #' 
  #' @description NATMI results can be easily thrown off if you don't use 
  #' specific file names. Overwriting, reading and appending to existing files 
  #' all at once can cause issues. Cross-talk caused from previous runs or 
  #' concurrent runs can create many errors and inaccuracies. Here we create 
  #' custom file paths for every input, output and intermediate value, to avoid
  #' this cross-talk. 
  #' 
  #' Using the systime and a user tag as unique markers, this function creates 
  #' NATMI-safe filepaths that can be directly passed to call_natmi.params.
  #' 
  #' @param tag An additional tag that specifies the save files beyond the time
  #' of run, as with parallel analyses the time can overlap.
  #' 
  #' @return A named list of custom filepaths.
  
  unique_natmi_filepaths <- function(tag) {
    
    # Time marker for this specific LIANA run. Underscores and hyphens can 
    # confuse natmi file recognition so we only use numbers.
    # We add a random five letter tag to further distinguish runs.
    natmi_output <-  Sys.time() %>%
      as.character()   %>%
      gsub(':', '', .) %>% 
      gsub(' ', '', .) %>%
      gsub('-', '', .) %>%
      str_glue('Job', tag, 'Test', .)
    
    # Use the unique tag to create unique filepaths
    natmi_params <-
      list(output_dir = natmi_output,
           expr_file  = str_glue(natmi_output, "_expr_matrix.csv"),
           meta_file  = str_glue(natmi_output, "_metadata.csv"),
           reso_name  = str_glue(natmi_output, "Resource"))
    
    
    return(natmi_params)
    
  }
}