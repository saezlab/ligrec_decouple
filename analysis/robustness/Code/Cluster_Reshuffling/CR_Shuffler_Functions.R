#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # This is the CR_Shuffler_Functions.R script. It defines all the core 
  # functions for reshuffling the cluster annotations of a given metadata. 
  # To see how it's implemented, check out the Run_Iterator.R script.
}



#------------------------------------------------------------------------------#
# 1. Define Functions ----------------------------------------------------------

# shuffle_Clusters()
{
  #' Reshuffle Cluster annotations in a Seurat
  #'
  #' @description This function takes a metadata table from a Seurat Object and
  #' partially reshuffles the cluster annotations of the cells in the table.
  #'
  #' N.B. that whenever an annotation is replaced, the replacement is drawn from
  #' the distribution of annotations already in the table that are not equal to
  #' the original annotation. This means that if n % of annotations are
  #' reshuffled, there will also be exactly n % mismatch between the old and new
  #' cluster annotations.
  #'
  #'
  #' @param seed As an integer. Cluster reshuffling has inherent
  #' randomness to it, by supplying a seed you ensure that this function runs
  #' reproducibly.
  #'
  #' @param mismatch_prop As a fraction. What proportion of the cluster
  #' annotations should be reshuffled? Also, what proportion of cluster
  #' annotations in the new metadata should differ from the old?
  #' (Both of these proportions are the same.)
  #'
  #' @param metadata The metadata df of a Seurat who's cell cluster annotations
  #' you want to reshuffle. Note that the metadata should have numeric cluster
  #' labels stored as a factor in a column called "cluster_key" of the metadata.
  #' The input metadata should also have the cell bar codes as its rownames.
  #'
  #'
  #' @return A new metadata data frame that has partially reshuffled/mismatched
  #' cluster annotations in the "cluster_key" column.
  
  
  shuffle_Clusters <- function(seed,
                               mismatch_prop,
                               metadata) {
    
    # There is randomness in reshuffling, so set the seed for reproducibility
    set.seed(seed)
    
    # We format the input df as a tibble and store that as metadata_old
    metadata_old <- metadata %>%
      rownames_to_column(var = "Bar_Code") %>%
      as_tibble()
    
    # Convert the cluster_key column froma  factor to an integer
    metadata_old$cluster_key <- metadata_old$cluster_key %>%
      as.numeric()
    
    
    
    
    # We select the rows we will reshuffle from the metadata. We keep the bar
    # code so we can remerge the old and reshuffled metadata later.
    meta_clusters_dilute <-
      slice_sample(metadata_old, prop = mismatch_prop) %>%
      select(c("Bar_Code", all_of("cluster_key")))
    
    # We map over the rows whose cluster_key assignment we want to reshuffle
    meta_clusters_dilute$cluster_key <-
      map(meta_clusters_dilute$cluster_key, function(annotation) {
        
        # We get the distribution of cluster annotations in metadata_old that 
        # are not equal to the annotation we are drawing from. 
        replacement_candidates <- metadata_old           %>%
          filter(metadata_old$cluster_key != annotation) %>%
          select(cluster_key)                            %>%
          unlist()
        
        # The new annotation is drawn at random from the above distribution. 
        # In this way we are guaranteed a mismatch through this replacement, but
        # in the big picture the proportion of clusters is not drastically 
        # altered.
        new_annotation <- sample(replacement_candidates, 1)
        
        return(new_annotation)
        
      }) %>%
      unlist() # because formatting
    
    
    
    # We join our reshuffled annotations with the old ones
    metadata_new <-
      left_join(metadata_old,
                meta_clusters_dilute,
                by = "Bar_Code",
                suffix = c("_old", ""))
    
    # This index marks every row that was not reshuffled
    helper_index <- which(is.na(metadata_new$cluster_key))
    
    # Carry over the original annotations from the unshuffled rows
    metadata_new$cluster_key[helper_index] <-
      metadata_new[["cluster_key_old"]][helper_index]
    
    # We add a column that marks which rows were reshuffled, and remove the old 
    # values from the tibble. Then we convert to a df because the input is a df.
    metadata_new <- metadata_new %>%
      mutate("Mismatched" = 
               if_else(.$cluster_key_old == .$cluster_key, FALSE, TRUE)) %>%
      select(-"cluster_key_old") %>%
      as.data.frame()
    
    
    
    # We convert back to the original format, renaming rownames, deleting the 
    # Bar_Code columns we created, converting the annotations back to factors
    rownames(metadata_new) <- metadata_new$Bar_Code
    
    metadata_new <- metadata_new %>%
      select(-"Bar_Code")
    
    metadata_new$cluster_key <- metadata_new$cluster_key %>%
      as.factor()
    
    
    
    # We give the user a status update
    print(str_glue("Cluster annotations reshuffled. ",
        round((sum(metadata_new$Mismatched) / nrow(metadata_new)) * 100, 2),
        " % mismatch to the original annotation."))
    
    
    
    # We return the results.
    return(metadata_new)
  
    
    } # end of function

  
}


# wrap_Shuffler()
{
  #' shuffle_Clusters Wrapper
  #' 
  #' @description This function prints a header to the console and then iterates
  #' the shuffle_Clusters() function over every seed provided for a given
  #' mismatch proportion.
  #' 
  #' @param seed_list A list of seeds to iterate over.
  #' 
  #' @param mismatch_prop A mismatch proportion for which to generate multiple
  #' permuations of reshuffled metadata.
  #' 
  #' @param metadata The metadata df of a Seurat who's cell cluster annotations
  #' you want to reshuffle. Note that the metadata should have numeric cluster
  #' labels stored as a factor in a column called "cluster_key" of the metadata.
  #' The input metadata should also have the cell bar codes as its rownames.
  #' 
  #' @return A list of reshuffled metadatas sorted by seed and mismatch.
  
  
  wrap_Shuffler <- function(seed_list,
                            mismatch_prop,
                            metadata) {
    
    # Let the user follow along in the console
    print_Title(str_glue(
      "Creating ",
      mismatch_prop * 100,
      " % mismatched cluster annotations."
    ))
    
    # Iterate the shuffler across all the seeds at this mismatch prop
    reshuffled_metadatas <- lapply(
      seed_list,
      shuffle_Clusters,
      mismatch_prop = mismatch_prop,
      metadata      = metadata
    )
    
    
    # Retrun our results.
    return(reshuffled_metadatas)
    
  } # end of function

  
  }
