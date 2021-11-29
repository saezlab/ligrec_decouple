#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  
  # This is the Dilutor_Processing_Functions.R script. These function work as 
  # intermediary steps in the process of diluting a resource. Information on
  # how exactly dilution is performed can be found in RD_Dilutor.R.
  
}



#------------------------------------------------------------------------------#
# 1. Defining Functions --------------------------------------------------------

# extract_unconflicting_Genes()
{
  #' Extract and filter genes
  #'
  #' @description Extracts a subset of genes from an input Seurat, and filters
  #' them to not include a given list of gene names.
  #'
  #' @param data_set The data set (as a Seurat object) from which you will draw
  #' genes from. Must be one using RNA assay data.
  #'
  #' @param feature_type Choose "generic" or "variable" (as a char).
  #'
  #' If generic, genes will be extracted from the row names of the RNA count
  #' matrix of the input, i.e. all the genes profiled in the data_set will be
  #' extracted.
  #'
  #' If variable extraction will be from the variable features of the data_set,
  #' which have to be set up beforehand.
  #'
  #' @param conflict_genes A list of genes that should not be in the output, as
  #' a list.
  #'
  #' @return Returns a list of genes from the data_set that does not include any
  #' members of the conflict_genes input.
  
  
  extract_unconflicting_Genes <- function (data_set,
                                           feature_type,
                                           conflict_genes) {
    # Depending on what type of dilution is requested, we pull our genes from
    # the variable or normal section of the seurat object.
    
    if (feature_type == "generic") {
      gene_name_list  <- as.list(rownames(data_set@assays$RNA@data))
      
      
    } else if (feature_type == "variable") {
      gene_name_list <- as.list(data_set@assays$RNA@var.features)
      
      
    } else {
      # Throw an error if feature_type is not "generic" or "variable".
      warning("FEATURE TYPE FOR DILUTION WAS NOT SET PROPERLY. RETURNING NULL.")
      return()
      
    }
    
    # Remove any blacklisted genes from gene_name_list
    gene_name_list <- gene_name_list %>%
      discard( ~ .x %in% conflict_genes)
    
    
    return(gene_name_list)
    
  } # end of function
  
  
}


# random_Dilute()
{
  #' Dilutes entire resources with a random method.
  #'
  #' @description Given a resource replaces all source_genesymbols and
  #' target_genesymbols with random relationships made from pairwise
  #' combinations of gene_name_list. In the process topological rules that mimic
  #' OmniPath are followed, namely:
  #'
  #' 1. No duplicate interactions created. All new interactions are unique.
  #' 2. In OmniPath, genes that are sources in some relationships and targets in
  #' others is rare. This is termed source and target overlap. All diluted
  #' interactions will have 0 % source and target overlap. A given gene will
  #' either always be a source, or always be a target.
  #' 3. No sources that are targets to themselves.
  #'
  #' Beyond these three rules interaction topology is random. This method
  #' requires less genes to work with than other methods, and if a solution in
  #' the above constraints is possible, this method will always find it.
  #' 
  #' The randomness within this function is controlled by the master_seed set in
  #' dilute_Resource().
  #'
  #' @param resource_dilute The resource (as a tibble) which you would like to
  #' falsify / dilute with random gene relationships. Must have a
  #' source_genesymbol and a target_genesymbol column. The method will introduce
  #' a column that marks all interactions in the diluted output as random. No
  #' other columns are modified in any way (LR_Pairs not updated for example).
  #'
  #' @param gene_name_list A list of gene names from which to generate random
  #' pairwise relationships.
  #'
  #' @return Returns a tibble that can be used as a custom OmniPath resource for
  #' liana_wrap and the call_method functions but has is entirely made up of
  #' random interactions, which are all marked.
  
  
  random_Dilute <- function (resource_dilute, gene_name_list) {
    
    # This method works by generating random interactions from gene_name_list.
    # Because OmniPath has a low overlap of sources and targets, i.e. few
    # gene products that are both ligands and receptor, we split
    # gene_name_list into sources and targets first. Since the highest number
    # of possible interactions in this scenario would come from equally sized
    # source and target lists, we split gene_name_list evenly.
    
    # Random sample for source_gene_list to avoid any ordering bias of
    # gene_name_list
    source_gene_list <- sample(gene_name_list,
                               size = ceiling(length(gene_name_list) / 2))
    
    # target_gene_list is what remains of gene_name_list when source_gene_list
    # is removed
    target_gene_list <- gene_name_list %>%
      discard( ~ .x %in% source_gene_list)
    
    
    
    
    # If we need to create more source/target interactions than is possible
    # given the gene_name_lists we have, we throw an error and don't even try
    # to compute it.
    max_possible_interactions <-
      length(source_gene_list) * length(target_gene_list)
    
    
    if (max_possible_interactions < nrow(resource_dilute)) {
      warning(
        str_glue(
          "THE NUMBER OF GENES GIVEN CANT CREATE AS MANY ",
          "RANDOM INTERACTIONS AS THE DILUTION PROP REQUIRES. ",
          "RETURNING NULL."
        )
      )
      return()
    }
    
    
    
    # We compute every possible unique interaction given our source and target
    # list and format the results. ST: Stands for Source and Target
    diluted_ST_interactions        <-
      expand.grid(source = source_gene_list,
                  target = target_gene_list)
    
    diluted_ST_interactions$source <-
      unlist(diluted_ST_interactions$source)
    diluted_ST_interactions$target <-
      unlist(diluted_ST_interactions$target)
    
    diluted_ST_interactions        <-
      tibble(diluted_ST_interactions)
    
    
    
    # We impute interactions into resource_dilute .
    
    # Because we are sampling from every possible unique interaction using the
    # setup where the maximum possible number of interactions is created, we
    # always find a solution if a solution with our criteria is possible.
    
    # We sample without replacing because we don't want duplicates.
    resource_dilute[c("source_genesymbol", "target_genesymbol")] <-
      slice_sample(diluted_ST_interactions, n = nrow(resource_dilute))
    
    
    # Mark every interaction as random, because every interaction is random.
    resource_dilute <- resource_dilute %>%
      mutate(isRandom = TRUE)
    
    # Return Output
    return(resource_dilute)
    
  } # end of function
  
  
}


# preserve_Dilute()
{
  #' Dilutes entire resources while semi-preserving topology
  #'
  #' @description Determines every unique gene in the given resource, then
  #' remaps every one of them on a 1:1 basis with genes provided as an argument.
  #' In the process, existing verified interactions from the resource are wholly
  #' replaced ("diluted") with false interactions. The resource is diluted. At
  #' the same time, since the process simply renames the genes, the resource
  #' topology is 100 % preserved. However, if the output was just a subset of
  #' a resource as in dilute_Resource(), the overall topology is only
  #' semi-preserved.
  #' 
  #' The randomness within this function is controlled by the master_seed set in
  #' dilute_Resource().
  #'
  #' @param resource_dilute The resource (as a tibble) which you would like to
  #' falsify / dilute with random gene relationships. Must have a
  #' source_genesymbol and a target_genesymbol column. The method will introduce
  #' a column that marks all interactions in the diluted output as random. No
  #' other columns are modified in any way (LR_Pairs not updated for example).
  #'
  #' @param gene_name_list A list of gene names to be used for dilution. The
  #' original gene names from resource_dilute will be wholly replaced with the
  #' gene names in this list.
  #'
  #' @return Returns a tibble that can be used as a custom OmniPath resource for
  #' liana_wrap and the call_method functions but has is entirely made up of
  #' diluted interactions, which are all marked.
  
  preserve_Dilute <- function (resource_dilute, gene_name_list) {
    
    # This method creates a dictionary, where every real gene in
    # resource_dilute gets a unique fake counterpart from gene_name_list.
    # Using the dictionary, we go through resource_dilute and swap each real
    # gene for its counterpart. As such, the resource_dilute topology is
    # preserved.
    
    # In order to create the dictionary, we extract all the unique gene names
    # in resource dilute
    resource_genes <-
      unlist(unique(
        c(
          resource_dilute$source_genesymbol,
          resource_dilute$target_genesymbol
        )
      ))
    
    # If gene_name_list isn't long enough to provide a 1:1 dictionary for
    # resource_genes we abort the process.
    if (length(gene_name_list) < length(resource_genes)) {
      warning(
        str_glue(
          "NOT ENOUGH GENES FROM DATA SET PROVIDED TO PRESERVE ",
          "TOPOLOGY. RETURNING NULL."
        )
      )
      return()
    }
    
    # We sample to avoid ordering bias in gene_name_list. We sample as many
    # genes as are in resource_genes to create our 1:1 dictionary.
    dilution_genes <- unlist(sample(gene_name_list,
                                    size = length(resource_genes)))
    
    # Define dictionary with resource_genes and dilution_genes
    dilution_dict <- tibble(resource_genes, dilution_genes)
    
    
    
    
    # Now that we have our dictionary, we iterate through every
    # source_genesymbol in resource_dilute, look up its counterpart in the
    # dictionary and replace source_genesymbol with the counterpart
    for (i in 1:nrow(resource_dilute)) {
      # Whats the name of the ource_genesymbol?
      gene_to_dilute  <- resource_dilute$source_genesymbol[i]
      
      # Where do I find it in the dictionary?
      GTD_index <-
        which(dilution_dict$resource_genes == gene_to_dilute)
      
      # Replace source_genesymbol with the counterpart at the row index
      # the source_genesymbol is at in the dictionary
      resource_dilute$source_genesymbol[i] <-
        dilution_dict[GTD_index, "dilution_genes"]
      
      
    }
    
    # Repeat the aboove process for target_genesymbol
    for (i in 1:nrow(resource_dilute)) {
      gene_to_dilute  <- resource_dilute$target_genesymbol[i]
      
      GTD_index <-
        which(dilution_dict$resource_genes == gene_to_dilute)
      
      resource_dilute$target_genesymbol[i] <-
        dilution_dict[GTD_index, "dilution_genes"]
      
      
    }
    
    # The above for loops nest resource_dilute, which we undo here. We also mark
    # all the interactions as random, because they are.
    resource_dilute <- resource_dilute %>%
      unnest(cols = c(source_genesymbol,  target_genesymbol)) %>%
      mutate(isRandom = TRUE)
    
    # Return Output
    return(resource_dilute)
    
  } # end of function
  
  
}


# dilution_Overview()
{
  #' Produces a console output for a diluted resource
  #' 
  #' @description  This function uses dilution parameters, the input resource
  #' and the results of dilution to output important statistics about how the 
  #' dilution process went to the user. While the statistics are output, the 
  #' warning logic also contains the information on whether any of the 
  #' displayed parameters are outside the bounds they are supposed to be in,
  #' alerting the user if something is wrong.
  #' 
  #' @param warning_logic A list derived from the dilute_Resource() function
  #' that states whether an error or warning in the dilution process has 
  #' occurred.
  #' 
  #' @param dilution_prop An integer indicating to what percent the input
  #' resource was supposed to be diluted to. 
  #' 
  #' @param feature_type Did dilution occur with all the genes profiled in
  #' testdata (choose "generic") or with the most variable features (choose 
  #' "variable")? As a string.
  #' 
  #' @param preserve_topology When diluting, two methods could've been
  #' used. random_Dilute makes only a small effort to preserve the 
  #' topology of the rows that are diluted from resource, preserve_Dilute makes 
  #' a far greater effort. Choose TRUE for preserve_Dilute and FALSE for 
  #' random_Dilute.
  #' 
  #' @param resource The original input resource before dilution. As a tibble.
  #' 
  #' @param new_resource The diluted resource derived from the original. As a 
  #' tibble.
  #' 
  #' @param resource_dilute The diluted rows that are in new_resource, as a 
  #' tibble.
  #' 
  #' @return Nothing directly, but an informative console output.
  
  
  dilution_Overview <- function(warning_logic,
                                
                                dilution_prop,
                                preserve_topology,
                                feature_type,
                                
                                resource,
                                new_resource,
                                resource_dilute) {
    

    # 1. Print Output Header (Recaps Dilution parameters)
    {
      # we briefly define this divider here so the code below is more readable
      output_divider <-
        "------------------------------------------------------------
      ------------------------------------------------------------"
      
      
      # spacing outputs helps with overview.
      print(str_glue(""))
      print(str_glue(""))
      
      # If an error occurred, this is the highest priority information. Thus 
      # this is what is brought to the users attention first. In this way, the 
      # output immediately warns the user of any problems and lets the user 
      # know which output is associated with which warnings.
      if (warning_logic$error_thrown == TRUE) {
        print(str_glue(
          as.character(dilution_prop * 100),
          "% Dilution --- ERROR OCCURRED"
        ))
        print(str_glue("Please check ERROR below."))
        
        
        # A warning is second priority. If there was an error that is more 
        # relevant than a warning, but if there was no error, a warning has the
        # next priority
      } else if (warning_logic$warning_thrown == TRUE) {
        print(str_glue(
          as.character(dilution_prop * 100),
          "% Dilution --- WARNING OCCURRED"
        ))
        print(str_glue("Please check WARNING below."))
        
        
        # Only if there were no errors or warnings do we label a dilution 
        # successful
      } else {
        print(str_glue(
          as.character(dilution_prop * 100),
          "% Dilution Successful"
        ))
        
      }
      
      print(str_glue(""))
      
      # The output reminds the user what the input parameters for dilution 
      # type were, and labels this output. When multiple outputs are produced
      # by iterating this function multiple times and the user is in doubt, 
      # they can hopefully piece together which dilution this is the output 
      # for.
      if (preserve_topology == TRUE) {
        # remind the user what "semi-preserved" means.
        print(
          str_glue(
            "Topology semi-preserved:
                       - Diluted interactions have identical topology.
                       - Genes NOT in dilution have identical topology.
                       - Genes that existed in both have split topology."
          )
        )
        
      } else if (preserve_topology == FALSE) {
        # remind the user of the topological constraints in random_Dilute()
        print(
          str_glue(
            "Topology unpreserved, but:
                       - No duplicate interactions.
                       - No sources that are their own target.
                       - No overlap between diluted sources and targets."
          )
        )
        
      }
      
      # space
      print(str_glue(""))
      
      # What did the feature type mean again?
      if (feature_type == "generic") {
        print(
          str_glue("Diluted with a list of all genes in the input Seurat."))
        
      } else if (feature_type == "variable") {
        print(str_glue(
          "Diluted with the variable features in the input ",
          "Seurat."
        ))
        
      }
    }
    
    
    
    
    # 2. Calculate lengthy expressions for output body:
    {
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
    }
    
    
    
    
    # 3. Print Output Body (Actual Resource Stats)
    {
      # Results for the analysis of the edge number and uniqueness
      print(str_glue(""))
      print(str_glue(output_divider))
      print(str_glue("Number of edges before:                           ",
                     nrow(resource)))
      print(str_glue("Number of edges after:                            ",
                     nrow(new_resource)))
      print(str_glue(""))
      
      print(str_glue("Number of unique edges before:                    ",
                     length(unique(
                       resource$LR_Pair
                     ))))
      print(str_glue("Number of unique edges after:                     ",
                     length(unique(
                       new_resource$LR_Pair
                     ))))
      print(str_glue(""))
      print(str_glue(output_divider))
      
      
      
      # Analysis of number of edges with LR_Pairs foreign to original resource
      # as well as the ones marked diluted. Also gives the achieves dilution
      # proportion compared the user defined one.
      print(str_glue(
        "Number of edges marked  as diluted (after):       ",
        sum(new_resource$isRandom)
      ))
      print(str_glue(""))
      
      print(str_glue(
        "Proportion of edges marked as diluted:            ",
        round(prop_marked_diluted, 2),
        " %"
      ))
      print(str_glue(
        "Proportion of edges actually diluted:             ",
        round(prop_marked_diluted, 2),
        " %"
      ))
      print(str_glue(""))
      print(str_glue(output_divider))
      
      
      
      # Limited analysis of topology in source + target overlap
      print(str_glue(
        "Source and Target overlap before:                 ",
        round(before_overlap, 2),
        " %"
      ))
      print(str_glue(
        "Source and Target overlap after:                  ",
        round(after_overlap, 2),
        " %"
      ))
      print(str_glue(""))
      print(str_glue(
        "Source and target overlap in diluted edges:       ",
        round(diluted_overlap, 2),
        " %"
      ))
      print(str_glue(""))
      print(str_glue(output_divider))
      
      
      
      # There should never be sources that are targets to themselves within a
      # single interaction.
      print(str_glue(
        "Sources that are Targets to themselves, before:   ",
        sum(
          resource$source_genesymbol ==
            resource$target_genesymbol
        )
      ))
      print(str_glue(
        "Sources that are Targets to themselves, after:    ",
        sum(
          new_resource$source_genesymbol ==
            new_resource$target_genesymbol
        )
      ))
      
      print(str_glue(""))
      print(str_glue(""))
      
    }
    
    
  } # end of function
    
  
}
