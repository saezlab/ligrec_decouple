#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # This is the Iterator_Top_Ranks.R script. It defines functions the Iterator 
  # uses to manipulate top-ranked CCIs.
}



#------------------------------------------------------------------------------#
# 1. Define Functions ----------------------------------------------------------

# clust_get_top_ranks()
{
  #' Get the top n ranked items of a method from the tibble liana wrapper or
  #' call_x results
  #'
  #' @param data_set The tibble output by the liana wrapper function or call_x
  #' functions.
  #' @param method The method for which you would like to extract the top ranked
  #' interactions, as a string.
  #' @param top_n The number of items to return, returns items ranked #1 to #n.
  #' @param with_ties TRUE or FALSE. An argument passed to slice_min and
  #' slice_max respectively. When slicing the top_n, sometimes a cutoff occurs
  #' between equally ranked items. Should the items also be included, expanding
  #' the selection (TRUE), or should preserving top_n be prioritized (FALSE).
  #'
  #' @return Returns the tibble input cut down to the top n highest ranked
  #' interactions.
  
  clust_get_top_ranks <-
    function(data_set, method, top_n, with_ties = FALSE) {
      # generate a list that describes how to rank in each method to get what the
      # method considers best
      rank_spec_list <-
        list(
          "cellchat"        = list(method_score = "pval",
                                   descending_order =  FALSE),
          
          "call_connectome" = list(method_score = "weight_sc",
                                   
                                   descending_order =  TRUE),
          
          "call_italk"      = list(method_score = "logfc_comb",
                                   descending_order =  TRUE),
          
          "call_natmi"      = list(method_score = "edge_specificity",
                                   descending_order =  TRUE),
          
          "call_sca"        = list(method_score = "LRscore",
                                   descending_order =  TRUE),
          
          "squidpy"         = list(method_score = "pvalue",
                                   descending_order =  FALSE)
        )
      
      
      
      # If its a p-value method, the code will go into this if statement
      if (rank_spec_list[[method]]$descending_order == FALSE) {
        topranks <-
          slice_min(
            data_set,
            n = top_n,
            order_by = !!sym(rank_spec_list[[method]]$method_score),
            with_ties = with_ties
          )
        
        # order_by is assigned the criterion dictated by rank_spec_list
        
        
      } else {
        # if it's not one of the p-value methods
        
        topranks <-
          slice_max(
            data_set,
            n = top_n,
            order_by = !!sym(rank_spec_list[[method]]$method_score),
            with_ties = with_ties
          )
        
        # order_by is assigned the criterion dictated by rank_spec_list
        
        
      }
      
      return(topranks)
      
    } #end of function
  
  
  
}


# format_top_ranks()
{
  #' Add an two columns, one contains unique Ligand-Receptor combinations, the 
  #' other unique Source-Target-Ligand-Receptor combinations.
  #' 
  #' @param top_ranks A tibble that has CCIs in it that have the source, target,
  #' ligand and receptor listed.
  #' 
  #' @return A tibble like the input with a column of unique tags per LR and
  #'  STLR.

  format_top_ranks <- function(top_ranks) {
    # Format the top_ranks data frame for future processing steps.
    top_ranks <- top_ranks %>%
      unite("LR_Pair",
            c(ligand, receptor),
            remove = FALSE,
            sep = "_") %>%
      relocate("LR_Pair", .after = last_col()) %>%
      unite("LR_ID",
            c(source, target, ligand, receptor),
            remove = FALSE,
            sep = "_") %>%
      relocate("LR_ID", .after = last_col())
    
    
    return(top_ranks)
    
  }
  

}



