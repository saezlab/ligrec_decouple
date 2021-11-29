# subset_Clusters()
{
  #' Subset Cluster annotations in a Seurat
  #'
  #' @description This function takes a metadata table from a Seurat Object and
  #' partially subsets the cluster annotations of the cells in the table.
  #'
  #'
  #' @param seed As an integer. Cluster subsetting has inherent randomness to 
  #' it, by supplying a seed you ensure that this function runs reproducibly.
  #'
  #' @param removal_prop A proportion of the cells to remove from each cluster.
  #' 
  #' @param metadata The metadata df of a Seurat who's cell cluster annotations
  #' you want to subset. Note that the metadata should have numeric cluster
  #' labels stored as a factor in a column called "cluster_key" of the metadata.
  #' The input metadata should also have the cell bar codes as its rownames.
  #'
  #'
  #' @return A new metadata data frame that has partially subsetted cluster 
  #' annotations in the "cluster_key" column.
  
  
  subset_Clusters <- function(seed,
                              removal_prop,
                              metadata) {
    
 
    # There is randomness in subsetting, so set the seed for reproducibility
    set.seed(seed)

    # We format the input df as a tibble and store that as metadata_old
    metadata_new <- metadata %>%
      rownames_to_column(var = "Bar_Code") %>%
      as_tibble() %>%
      group_by(cluster_key) %>%
      slice_sample(prop = 1-removal_prop) %>%
      ungroup() %>%
      as.data.frame()
    
    # We convert back to the original format, renaming rownames, deleting the
    # Bar_Code columns we created, converting the annotations back to factors
    rownames(metadata_new) <- metadata_new$Bar_Code

    metadata_new <- metadata_new %>%
      select(-"Bar_Code")

    metadata_new$cluster_key <- metadata_new$cluster_key %>%
      as.factor()
    
    
    # Calculate cluster proportions difference and siplay to user.
    {
      
      cluster_counts <- metadata %>%
        count(cluster_key) %>%
        rename("old_counts" = n) %>%
        mutate("new_counts" = count(metadata_new, cluster_key)$n)
      
      cluster_props <- cluster_counts %>%
        mutate("old_proportions" = old_counts / sum(old_counts)) %>%
        mutate("new_proportions" = new_counts / sum(new_counts)) %>%
        mutate("difference" = abs(old_proportions - new_proportions)) %>%
        select(-c("old_counts", "new_counts")) 
      
      cluster_props[, -1] <- round(cluster_props[, -1], 3)
      
      cluster_props <- cluster_props %>% as_tibble()
      
      
      # We give the user a status update
      cat("Cluster annotations were subset. Proportions comparison: \n\n")
      print(cluster_props)
      cat("\n\n\n")
      
      
    }


    
    # We return the results.
    return(metadata_new)
    
    
  } # end of function
  
  
}


# wrap_Subsetter()
{
  #' subset_Clusters Wrapper
  #' 
  #' @description This function prints a header to the console and then iterates
  #' the subset_Clusters() function over every seed provided for a given
  #' mismatch proportion.
  #' 
  #' @param seed_list A list of seeds to iterate over.
  #' 
  #' @param removal_prop A proportion of the cells to remove from each cluster.
  #' 
  #' @param metadata The metadata df of a Seurat who's cell cluster annotations
  #' you want to subset. Note that the metadata should have numeric cluster
  #' labels stored as a factor in a column called "cluster_key" of the metadata.
  #' The input metadata should also have the cell bar codes as its rownames.
  #' 
  #' @return A list of subset metadatas sorted by seed and mismatch.
  
  wrap_Subsetter <- function(seed_list,
                             removal_prop,
                             metadata) {
    
    # Let the user follow along in the console
    print_Title(str_glue(
      "Creating ",
      removal_prop * 100,
      " % subset cluster annotations."
    ))
    
    # Iterate the Subsetter across all the seeds at this mismatch prop
    subset_metadatas <- lapply(
      seed_list,
      subset_Clusters,
      removal_prop  = removal_prop,
      metadata      = metadata
    )
    
    
    # Retrun our results.
    return(subset_metadatas)
    
  } # end of function
  
  
}