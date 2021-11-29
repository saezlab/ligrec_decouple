#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  
  # This is the Iterator_Processing_Functions.R script. These function work as 
  # miscellaneous intermediary steps in the process of running the iterator. 
  # Information on how the iterator runs can be found in 
  # RD_Robustness_Iterator.R.
  
}



#------------------------------------------------------------------------------#
# 1. Defining Functions --------------------------------------------------------

# reformat_Results()
{
  #' Reformats the outputs from lapplying resource_Robustness over the master
  #' seed list into a better organized hierarchy
  #' 
  #' @description The standard outputs from lapplying resource_Robustness() are 
  #' highly nested and unintuitive, this being a product of how they are 
  #' created. This function uses transpositions, flattenings and renaming to 
  #' create a flatter and easier to understand result object from iterating 
  #' resource_Robustness(), which can be used for further data extraction more 
  #' easily.
  #' 
  #' @param results The output from lapplying resource_Robustness over the 
  #' master seed list.
  #' 
  #' @return A better structured version of the input.
  
  reformat_Results <- function(results) {
    
    # At most, there are six outputs. They fall into three pairs of two that are
    # each formatted the same way. Liana results and top ranks are formatted
    # the same way, resources and top_ranks analysis are formatted the same way,
    # and runtime and testdata is formatted the same way.
    
    
    # We start by transposing results
    results <- transpose(results)
    
    
    # This is the segment of the results containing runtime and testdata
    segment_runtime_test <-
      results[names(results) %in% intersect(names(results),
                                            c("runtime",
                                              "testdata"))]
    
    # This is the segment of the results containing resources and analysis
    segment_resources_analysis <-
      results[names(results) %in% intersect(names(results),
                                            c("resources_OP",
                                              "top_ranks_analysis"))]
    
    # This is the segment of the results containing ranks and results
    segment_results_ranks <-
      results[names(results) %in% intersect(names(results),
                                            c("liana_results_OP",
                                              "top_ranks_OP"))]
    
    # Runtime_test is already corrrectly formatted.
    
    # Format resources_analysis, we transpose at a deeper level
    segment_resources_analysis <- segment_resources_analysis %>%
      map_depth(.depth = 1, transpose) %>%
      flatten_names(depth = 1)
    
    # Format results_ranks
    segment_results_ranks <- segment_results_ranks %>%
      map_depth(.depth = 1, transpose) %>%
      map_depth(.depth = 2, transpose) %>%
      flatten_names(depth = 2)
    
    # Recombine our results again and return them.
    restructured_results <- list(segment_results_ranks, 
                                 segment_resources_analysis,
                                 segment_runtime_test) %>%
      flatten()
    
    return(restructured_results) 
    
  }
}


# flatten_names()
{
  #' A helper function for restructure results, a variant on flatten() that 
  #' preserves names
  #' 
  #' @description This function takes a list with minimum three layers (a 
  #' sublist and a sub-sub list) and renames the lowest elements to include the 
  #' name of their parent list one step up in the hierarchy, and then flattens 
  #' that level of the list. In short, it removes one level of hierarchy from 
  #' the bottom, but makes sure that the names of the items indicate what 
  #' sublist they used to be a part of. 
  #' 
  #' This helps avoid situations where the lowest tier items in a list all have 
  #' the same names, the list gets flattened and you can't tell where the items
  #' came from anymore. 
  #' 
  #' @param high_tier_list The list to flatten the ends of.
  #' 
  #' @param depth The depth of the lowest items in the list. 
  #' 
  #' @return The input list with it's lowest level of hierarchy removed, and the
  #' lowest level items in it renamed to show whqat sublist they used to be in.
  
  
  flatten_names <- function(high_tier_list, depth) {
    
    new_high_tier <-
      # map to the user specififed depth
      map_depth(high_tier_list, depth, function(two_tier_list) {
        
        # Once at the second lowest tier, grab the name of the lists here
        # and the lists themselves.
        # We then map to the lowest tier, with the names in hand and rename the 
        # lowest tier elements.
        two_tier_list %>%
          map2(names(.), function(one_tier_list, one_tier_list_name)
            rename_list(one_tier_list, one_tier_list_name)) %>%
          flatten()
        
      })
    
    # return our new list
    return(new_high_tier)
    
  }
}


# rename_list()
{
  #' Helper function for reformat_Results, takes a list element and adds a tag 
  #' to its name
  #' 
  #' @param list_element A list element to rename.
  #' 
  #' @param str The tag you want to slap onto the list elements name.
  #' 
  #' @description Usually this is passed the lowest elements in a nested list
  #' as list_elements and then the name of the sub-list they are a part of as
  #' str.
  #' 
  #' @return A renamed list_element.
  
  rename_list <- function(list_element, str){
    
    new_list <- setNames(list_element,
                         str_glue("{str}_{names(list_element)}"))
    
    return(new_list)
    
  }
}


# extract_top_ranks()
{
  #' Extracts top_rank_overlaps from the restructured resource_Robustness() 
  #' output and formats and collates it into one convenient tibble
  #' 
  #' @param results The restructured resource_Robustness() output, i.e. the 
  #' output from reformat_Results().
  #' 
  #' @return A tibble that collates all the top_ranks overlap data from the 
  #' input.
  
  extract_top_ranks <- function(results) {
    
    # get just the top_ranks anylsis information
    top_ranks_analysis <- results$top_ranks_analysis
    
    # narrow it down to top_ranks_overlaps tibbles by searching for entries that 
    # contain the word "overlap"
    where_overlap <- str_detect(names(top_ranks_analysis), "Overlap")
    
    # bind all these separate tibbles into one 
    collated_top_ranks_overlap <- top_ranks_analysis[where_overlap] %>%
      bind_rows() 
    
    # reorganize the tibble to be more usable and easily understandable 
    collated_top_ranks_overlap <- collated_top_ranks_overlap %>%
      arrange(dilution_prop) %>%
      pivot_longer(cols = !(starts_with("dilution_prop")), names_to = "Method") %>%
      arrange(Method) %>%
      rename("Overlap" = value) 
    
    # remove unnecessary clutter
    rm(where_overlap)
    
    # returnt eh collated top_ranks_overlap information
    return(collated_top_ranks_overlap)
  }
}


# auto_plot_Description()
{
  #' Automatically creates a verbose caption of a top ranks overlap plot
  #' 
  #' @param top_ranks_overlap As a tibble in the form of a extract_top_ranks 
  #' output, though ideally it will be preprocessed for plotting (better method 
  #' names, no NAs, etc.). This is the top_ranks_overlap that would be plotted
  #' with this caption. The function takes data from the tibble's general 
  #' structure to describe it accurately. 
  #' 
  #' @param trial_run The same parameter from wrap_resource_Iterator(). Used
  #' in the file name to mark the file.
  #' 
  #' @param modify_baseline TRUE or FALSE. Should the top-ranked interactions from
  #' the baseline be modifiable by dilution. Usually this is not the case, and 
  #' almost all documentation is written from the perspective that the baseline is
  #' not modifiable.
  #'
  #' @param preserve_topology The same parameter from 
  #' wrap_resource_Iterator(). Used in the file name to mark the file.
  #' 
  #' @param testdata_type The same parameter from wrap_resource_Iterator(). 
  #' Used in the file name to mark the file.
  #' 
  #' @param feature_type The same parameter from wrap_resource_Iterator(). 
  #' Used in the file name to mark the file.
  #' 
  #' @param number_ranks The same parameter from wrap_resource_Iterator(). 
  #' Used in the file name to mark the file.
  #' 
  #' @param time_of_run The char tag of the time the script started being 
  #' executed.
  #' 
  #' @return A verbose caption describing the parameters used to generate the 
  #' results in the plot.
  
  auto_plot_Description <- function(top_ranks_overlap,
                                    
                                    trial_run,
                                    preserve_topology,
                                    modify_baseline,
                                    testdata_type,
                                    feature_type,
                                    number_ranks,
                                    time_of_run) {
    
    ## General comment, on testdata type, feature_type and topology, and if
    ## the baseline was modifiable
    {
      if (preserve_topology == FALSE) {
        topology_comment <- "random_Dilute()"
        
      } else if (preserve_topology == TRUE) {
        topology_comment <- "preserve_Dilute()"
        
      }
      
      
      if (modify_baseline == FALSE) {
        baseline_comment <- 
          "The baseline CCC predictions were preserved in the resource."
        
      } else if (modify_baseline == TRUE) {
        baseline_comment <- 
          "The baseline CCC predictions were modified in the resource. "
        
      }
      
      
      general_comment <-
        str_glue(
          "This plot was created using the ",
          testdata_type,
          " data. Dilution was performed using ",
          feature_type,
          " features and the ",
          topology_comment,
          " function. \n",
          baseline_comment
        )
      
      rm(topology_comment)
    }
    
    
    
    ## Dilution comment, on proportions
    {
      dilution_overview <- count(top_ranks_overlap,
                                 dilution_prop,
                                 run_mode = "real")
      
      
      dilution_comment <- str_glue(
        "The dilution occured in ",
        dilution_overview$dilution_prop[2] -
          dilution_overview$dilution_prop[1],
        " % increments up to a maximum of ",
        max(top_ranks_overlap$dilution_prop),
        " %. "
      )
      
      if (nrow(dilution_overview) < 1) {
        stop(
          "Expected at least two dilution proportions in input (0, and one ",
          "more. But found only one instead, namely ",
          dilution_overview$dilution_prop
        )
      }
      
      if (length(unique(dilution_overview$n)) != 1) {
        stop(
          "There should be an equal number of samples for every dilution, ",
          "but there is not."
        )
      }
      
      
      rm(dilution_overview)
      
    }
    
    
    ## Nperms and top_ranks comment
    {
      top_ranks_vector <- unlist(number_ranks)
      
      permutations_overview <- top_ranks_overlap %>%
        filter(dilution_prop == 0) %>%
        count(Method)
      
      
      top_ranks_permutations_comment <-
        str_glue(
          "The overlap was compared between the ",
          median(top_ranks_vector),
          " highest ranked interactions over ",
          permutations_overview$n[1],
          " permutations."
        )
      
      if (length(unique(permutations_overview$n)) != 1) {
        stop(
          "There should be an equal number of samples for each method at , ",
          "dilution proportion 0, but there is not."
        )
      }
      
      rm(permutations_overview, top_ranks_vector)
    }
    
    
    ## Date and time comment
    time_comment <- str_glue("Generated at ",
                             time_of_run,
                             ".")
    
    
    ## Assemple plotting caption
    plotting_caption <-
      str_glue(
        general_comment,
        "\n",
        dilution_comment,
        "\n\n",
        top_ranks_permutations_comment,
        "\n",
        time_comment
      )
    
    
    ## Add addendum if trial run
    if (trial_run == TRUE) {
      plotting_caption <-
        str_glue(plotting_caption, "   --   [TRIAL RUN]")
    }
    
    return(plotting_caption)
  }  # end of function
  
}

