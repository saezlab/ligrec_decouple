#------------------------------------------------------------------------------#
# 0. Define dilute_Resource() --------------------------------------------------

# dilute_Resource()
{
  #' Dilutes a subset of a resource while analysing its output
  #'
  #' @description This function acts as a handler for the entire process of
  #' diluting a resource. Diluting a resource is the process of replacing
  #' interactions of a resource with new interactions not from the resource.
  #' These new interactions are obviously not supported by the literature
  #' the resource is based on, they represent spurious "fake" interactions.
  #' From the viewpoint of the resource, they are unverified or irrelevant to
  #' CCI.
  #'
  #' Adding these interactions can model the addition of lower quality noise
  #' to a resource or mimic the switching from one resource to another. A
  #' certain core remains the same while the rest is exchanged. When we rerun
  #' liana methods with a diluted resource, the degree of similarity of its
  #' predictions gives insight into the robustness of the method as the
  #' resource changes. In particular,gradual changes in resource can be
  #' achieved that aren't possible when a resource is swapped wholesale.
  #'
  #' This function first formats an input resource so that only a proportion
  #' (dilution_prop) of it is diluted. Optionally, a list of interactions
  #' (top_rank_list) can be supplied, when the resource is diluted these
  #' interactions will be untouched. Commonly this is used so that top_ranked
  #' interactions aren't altered by dilution, so that only the false positive
  #' rate is measured. Ideally, the input resource should only include
  #' interactions that are also possible in data_set. Excess rows only
  #' increase the chances that dilution un-uniformly affects the rows that are
  #' relevant.
  #'
  #' Once the subset of the resource that will be diluted is known, a list of
  #' genes to dilute with is necessary. Ideally , these genes should come from
  #' the data set you will later use in LIANA with the diluted resource you
  #' generate. Since the resource only contains relevant genes, it is only
  #' accurate to replace them with genes that will also be relevant. To do this,
  #' the method runs  extract_unconflicting_Genes() on data_set. It will extract
  #' genes of the given feature_type.
  #'
  #' With these prerequesites, the method calls either random_Dilute() or
  #' preserve_Dilut() to perform the actual dilution of the subset using the
  #' list of genes from extract_unconflicting_Genes().
  #'
  #' Finally, the method analyses the diluted output resource in detail and
  #' outputs warnings if anything is off about it. A detailed output can also
  #' be achieved with the verbose argument.
  #'
  #' @param resource The resource (as a tibble) which you would like to
  #' falsify / dilute with new gene relationships. Ideally, the resource
  #' should only include interactions that are also possible in the data set.
  #' Dilution should uniformly affect the relevant parts of a resource after
  #' all, and not the rest.
  #'
  #' @param top_rank_list Optional. Gene interactions that should not be
  #' diluted, given as a list of chars or vector of chars. Every item of the
  #' list or vector should be formatted as an LR_Pair (source + _ + target,
  #' e.g "ITGB2_ICAM2"). Often this will be the LR_Pair column of a top_ranks
  #' data frame. For example, if the output of get_top_n_ranks() was
  #' top_ranks_df, top_rank_list = top_ranks_df$LR_Pair would work.
  #'
  #' @param dilution_prop As a number between 0-1. The proportion of rows of the
  #' resource to dilute. Interactions in top_rank_list can't be diluted. If
  #' attaining the requested dilution proportion requires overwriting top ranked
  #' interactions, the function throws an error and returns nothing instead.
  #'
  #' @param preserve_topology Either TRUE or FALSE. If FALSE, uses
  #' random_Dilute() to dilute the resource. If TRUE,uses preserve_Dilute() to
  #' dilute the resource.
  #'
  #' @param data_set A parameter for extract_unconflicting_Genes().
  #' Specifically, the data_set the gene names for dilution are extracted from.
  #'
  #' @param feature_type A parameter for extract_unconflicting_Genes(). More
  #' specifically, the type of genes to extract from data_set. More information
  #' can be found in the extract_unconflicting_Genes() documentation.
  #'
  #' @param verbose Set to TRUE by default. Produces a detailed output that
  #' summarizes important elements of the resource and how they changed through
  #' dilution. The output can be used to double check the parameters that were
  #' used and whether or not the dilution went correctly.
  #'
  #' Normally, the number of unique edges should not change, all edges should be
  #' unique, the percentage of edges MARKED diluted should be equal to the
  #' percentage that ARE diluted and should be close to the given dilution
  #' proportion. The overlap in diluted edges should be 0 if
  #' preserve_topology == FALSE, otherwise it should be close to the source and
  #' target overlap before dilution. No sources should be targets to themselves,
  #' ever.
  #'
  #' @param master_seed  At some stages when diluting  a resource, randomness
  #' is at play. Setting a master seed will ensure that the function runs the
  #' same way every time. Each instance of randomness within the scope of this
  #' function will run reproducibly.
  #'
  #' @return Returns a tibble that can be used as a custom OmniPath resource for
  #' liana_wrap and the call_method functions but has a certain (marked)
  #' percentage of it diluted with  interactions.
  
  
  dilute_Resource <- function(resource,
                              top_rank_list = list(),
                              dilution_prop,
                              preserve_topology,
                              data_set,
                              feature_type,
                              verbose = TRUE, 
                              master_seed)         {
    
    ## 0. Setting the Seed
    {
      # We set the seed once here, all functions using randomness within the scope
      # of this function will now run reproducibly.
      set.seed(master_seed)
    }
    
    
    ## 1. Splitting the Resource
    {
      # We divide the resource into three separate subsections.
      #   1. Interactions that are in top_rank_list, i.e. ones that should not 
      #      be diluted, named resource_top
      #   2. Interactions that are free to be diluted, named resource_bottom
      #   3. The specific interactions from resource_bottom that will be diluted
      #      to most closely match dilution_prop
      
      
      # Separate top-ranked parts of resource, create resource_top
      resource_top <- resource %>%
        filter(LR_Pair %in% top_rank_list)
      
      # Everything that is not top_ranked is open for dilution
      resource_bottom <- resource %>%
        filter(!(LR_Pair %in% top_rank_list))
      
      # Determine how many rows of resource_bottom to dilute so that the overall
      # dilution_prop is met.
      dilution_number <- round(nrow(resource) * dilution_prop)
      
      # Additional Warning message and break if the dilution prop can't be met
      # without diluting top-ranked interactions (which this function doesn't do)
      if (dilution_number > nrow(resource_bottom)) {
        warning(
          str_glue(
            "REQUESTED DILUTION PROPORTION NOT ATTAINABLE WITHOUT ",
            "OVERWRITING THE GIVEN TOP-RANKED INTERACIONS. ",
            "RETURNING NOTHING INSTEAD."
          )
        )
        
        return()
        
        
      }
      
      
      # If we dilute dilution_number of rows dilution_prop will be met
      # Copy dilution_number rows at random from resource_bottom...
      resource_dilute <-
        slice_sample(resource_bottom, n = dilution_number)
      
      # ... and delete them from resource bottom.
      resource_bottom <- anti_join(resource_bottom, resource_dilute)
      
      # Currently, resource_top + resource_bottom + resource_dilute == resource
      # Once resource_dilute has been properly diluted with the methods below
      # resource_top + resource_bottom + resource_dilute* == new_resource
      # The number of rows of new_resource that are diluted are equal to
      # dilution_prop.
      
    } # end of subpoint
    
    
    ## 2. Preparing a Gene Name List
    {
      # As described above, we need a gene_name_list to create diluted
      # interactions with. The extract_unconflicting_Genes() function extracts
      # genes from the data_set and removes any overlap to existing genes in the
      # resource.
      
      gene_name_list <-
        extract_unconflicting_Genes(
          data_set       = data_set,
          feature_type   = feature_type,
          conflict_genes =
            c(
              resource$source_genesymbol,
              resource$target_genesymbol
            )
        )
      
      # With the dilution candidates in resource_dilute and dilution genes 
      # ready, we move on to the actual dilution methods.
      
    } # end of subpoint
    
    
    ## 3. Implementing Dilution Method
    {
      # There are two approaches to dilution. One is random and that follows
      # certain topological rules but makes no further effort to match the
      # topology of OmniPath (yes, these rules are specific to OmniPath and not
      # generalization to other resources). The other preserves the topology of
      # resource_dilute exactly, and thus semi-preserves the topology of 
      # resource in new_resource. More details in the documentation.
      
      # random_Dilute(). There will be no duplicates, no sources that are their
      # own target, and no source and target overlap.
      if (preserve_topology == FALSE) {
        resource_dilute <- random_Dilute(
          resource_dilute = resource_dilute,
          gene_name_list  = gene_name_list
        )
        
      }
      
      # preserve_Dilute(). The topology of resource_dilute will be unaltered.
      if (preserve_topology == TRUE)  {
        resource_dilute <-
          preserve_Dilute(
            resource_dilute = resource_dilute,
            gene_name_list  = gene_name_list
          )
        
      }
      
      
      # Update the LR-Pair column with the new random "interaction partners".
      resource_dilute <- resource_dilute %>%
        select(-LR_Pair)                 %>%
        unite("LR_Pair",
              c(source_genesymbol, target_genesymbol),
              remove = FALSE)            %>%
        relocate("LR_Pair", .after = last_col())
      
      
      # The new resource has top ranked interactions, non-top rank but still 
      # real interactions, and diluted random interactions.
      new_resource <-
        bind_rows(resource_top, resource_bottom, resource_dilute)
      
    } # end of subpoint
    
    
    ## 4. Output Analysis and Warnings Assessment
    {
      # At this point, the output is ready, but to make sure nothing has gone
      # wrong in the dilution process, we calculate multiple metrics to capture
      # the topology of the result.
      
      # We start by calculating some more complicated expressions to save
      # typespace elsewhere.
      
      # calculate what proportion of edges is marked as diluted
      prop_marked_diluted <-
        sum(new_resource$isRandom) * 100 / nrow(new_resource)
      
      # calculate what proportion of edges have an LR_Pair unique to the new
      # resource
      prop_diluted <- sum(!(new_resource$LR_Pair %in%
                              resource$LR_Pair))   * 100 / nrow(new_resource)
      
      # calculate percentage of source & target overlap in the old resource,
      # new resource, and in the specific rows that were diluted.
      before_overlap <-
        sum(resource$source_genesymbol %in%
              resource$target_genesymbol)          * 100 / nrow(resource)
      
      after_overlap <-
        sum(new_resource$source_genesymbol %in%
              new_resource$target_genesymbol)      * 100 / nrow(new_resource)
      
      diluted_overlap <-
        sum(resource_dilute$source_genesymbol %in%
              resource_dilute$target_genesymbol)   * 100 / nrow(resource_dilute)
      
      
      # Having our complex parameters sorted, we now start a series of 
      # evaluations as to whether or not the topology looks healthy. Guidelines 
      # for this are also in the documentation of the verbose argument above.
      
      # Store our evlauation results in a compact format
      warning_logic <- list()
      
      # The number of edges should be constant, and every edge should be unique
      warning_logic["edges_okay"]       <-
        all(sapply(list(
          nrow(resource),
          nrow(new_resource),
          length(unique(resource$LR_Pair)),
          length(unique(new_resource$LR_Pair))
        ),
        function(x)
          x == nrow(resource)))
      
      # The number of rows marked with isRandom == TRUE should be equal to the
      # number of rows with LR_Pairs foreign to the original resource
      warning_logic["dil_marker_okay"]  <-
        prop_marked_diluted == prop_diluted
      
      # The proportion of rows diluted shouldn't be too different from the user
      # defined proportion
      warning_logic["dil_prop_okay"]    <-
        abs(dilution_prop * 100 - prop_diluted) < 5
      
      # If random_Dilute() was used, the source and target overlap should be 0
      warning_logic["r_topology_okay"]  <-
        diluted_overlap == 0
      
      # If preserve_Dilute() was used, the source and target overlap should be
      # similar to the original input
      warning_logic["p__topology_okay"] <-
        abs(before_overlap - after_overlap) <= 10
      
      # There should never be sources that are targets to themselves in an
      # interaction
      warning_logic["no_self_targets"]  <-
        all(sapply(list(
          sum(
            resource$source_genesymbol ==
              resource$target_genesymbol
          ),
          sum(
            new_resource$source_genesymbol ==
              new_resource$target_genesymbol
          )
        ),
        function(x)
          x == 0))
      
      
      
      
      # In this next if statement we integrate the above knowledge, is there an
      # issue in the data? If so, is it mild (a warning) or severe (an error).
      
      # We only want to consider r_topology_okay if random_Dilute() was used
      if (preserve_topology == FALSE) {
        # The dilution proportion being off is concerning but not necessarily
        # indicative that the dilution went wrong. We want to distinguish this
        # issue as a warning, not an error.
        warning_logic["warning_thrown"]   <-
          !(all(warning_logic$dil_prop_okay))
        
        # Any if these issues are serious and indicate the dilution failed.
        warning_logic["error_thrown"]     <-
          !(
            all(
              warning_logic$edges_okay,
              warning_logic$dil_marker_okay,
              warning_logic$r_topology_okay,
              warning_logic$no_self_targets
            )
          )
        
        # We only want to consider p_topology_okay if preserve_Dilute() was used
      } else if (preserve_topology == TRUE) {
        # As above, these two issues are concerning when they occur, but only
        # at the level of a warning. The dilution may still be fine.
        warning_logic["warning_thrown"]   <-
          !(all(
            warning_logic$dil_prop_okay,
            warning_logic$p_topology_okay
          ))
        
        # As above, these issues are serious and indicate the dilution failed.
        warning_logic["error_thrown"]     <-
          !(
            all(
              warning_logic$edges_okay,
              warning_logic$dil_marker_okay,
              warning_logic$no_self_targets
            )
          )
        
        # remove the cluttered parameters now that warnings are assessed warnings
        rm(
          before_overlap,
          after_overlap,
          diluted_overlap,
          prop_marked_diluted,
          prop_diluted,
          output_divider
        )
        
      }
      
      
      
    } # end of subpoint
    
    
    ## 5. Optional Verbose Output
    {
      # In this segment we print a bunch of information on key parameters of the
      # the inputs and outputs using the dilution_Overview() function.
      
      if (verbose == TRUE) {
        
        dilution_Overview(warning_logic     = warning_logic,
                          
                          dilution_prop     = dilution_prop,
                          preserve_topology = preserve_topology,
                          feature_type      = feature_type,
                          
                          resource          = resource,
                          new_resource      = new_resource,
                          resource_dilute   = resource_dilute)
        
        
      }
      

      
    } # end of subpoint
    
    
    ## 6. Print Warnings
    {
      # We print the warnings here even though we've had all the knowledge on
      # which warnings are needed for a while. We print them here because it is
      # more intuitive for the warnings to be right under the verbose output.
      
      if (warning_logic$edges_okay == FALSE) {
        warning(
          str_glue(
            "DILUTION ERROR: THE INPPUT AND OUTPUT SHOULD ONLY ",
            "INCLUDE UNIQUE INTERACTIONS AND SHOULD HAVE THE ",
            "SAME NUMBER OF INTERACTIONS. THIS IS NOT THE CASE."
          )
        )
        
      }
      
      
      if (warning_logic$dil_marker_okay == FALSE) {
        warning(
          str_glue(
            "DILUTION ERROR: THE PROPORTION OF INTERACTIONS ",
            "MARKED AS DILUTED (isRandom) DOES NOT CORRESPOND TO ",
            "THE PROPORTION OF INTERACTIONS THAT ACTUALLY ARE DILUTED. ",
            "SOMETHING IS WRONG WITH THE isRandom MARKINGS."
          )
        )
      }
      
      
      if (warning_logic$dil_prop_okay == FALSE) {
        warning(
          str_glue(
            "DILUTION WARNING: THE ACTUAL DILUTION PROPORTION IS ",
            "MORE THAN 5 PERCENTAGE POINTS DIFFERENT TO THE USER ",
            "SPECIFIED DILUTION PROPORTION. SOMETHING IS WORNG ",
            "WITH THE NUMBER OF DILUTED ROWS."
          )
        )
      }
      
      
      # Only consider r_topology_okay if random_Dilute() was used
      if (preserve_topology == FALSE) {
        if (warning_logic$r_topology_okay == FALSE) {
          warning(
            str_glue(
              "DILUTION ERROR: THE SOURCE AND TARGET  OVERLAP ",
              "SHOULD BE 0 FOR THIS DILUTION TYPE AND ISN'T. THE ",
              "TOPOLOGY IS NOT AS IT SHOULD BE."
            )
          )
        }
        
        
        # Only consider p_topology_okay if preserve_Dilute() was used
      } else if (preserve_topology == TRUE) {
        if (warning_logic$p__topology_okay == FALSE) {
          warning(
            str_glue(
              "DILUTION WARNING: THE SOURCE AND TARGET  OVERLAP ",
              "CHANGED MORE THAN 10 PERCENTAGE POINTS THROUGH ",
              "DILUTION. THE TOPOLOGY HAS CHANGED MORE ",
              "DRASTICALLY THAN EXPECTED FOR THIS METHOD."
            )
          )
        }
        
        
      }
      
      
      if (warning_logic$no_self_targets == FALSE) {
        warning(
          str_glue(
            "DILUTION ERROR: THE ORIGINAL RESOURCE AND/OR THE ",
            "DILUTED RESOURCE CONTAINS INTERACTIONS WITH ",
            "SOURCES THAT ARE TARGETS TO THEMSELVES."
          )
        )
      }
      
      rm(warning_logic)
      
      
    } # end of subpoint
    
    
    ## 7. Return Output
    return(new_resource)
    
    
  } #end of function
  
  
  
}


