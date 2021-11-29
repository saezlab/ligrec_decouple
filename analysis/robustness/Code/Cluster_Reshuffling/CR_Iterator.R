#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # This is the CR_Iterator.R script. It serves as the core structure that runs
  # the robustness analysis mentioned in Run_Iterator.R. It does this by 
  # defining the core wrap_cluster_Iterator() function, which goes through all
  # the necessary analysis steps. 
  # 
  # To get an overview of the conceptual approach, check the Run_Iterator.R 
  # script.
  # 
  # To get an overview of the programmatic approach, read through the 
  # subsections of the function definition below. The Iterator function also 
  # makes use of many sub-functions, whose documentation and code can be found
  # in other scripts from the code folder.
  
}



#------------------------------------------------------------------------------#
# 1. Define wrap_cluster_Iterator() -----------------------------------------

#' Determine the robustness of CCI inference methods in regards to cluster 
#' annotation accuracy
#'
#' @description This function evaluates the robustness of CCI-Inference methods'
#' predictions. It does so by running them on default cluster annotations and 
#' reshuffled cluster annotations and then comparing the overlap between the 
#' methods. This shows the methods ability to recreate the original input in
#' relation to the degree of mismatch between the original and reshuffled 
#' annotations. Besides a full analysis, this function also records metadata,
#' plots results and optionally saves results to the outputs folder.
#' 
#' @param testdata A Seurat containing scRNA data that has been preprocessed.
#' 
#' @param testdata_type A string that serves as a name or label for your 
#' testdata. It will be included in the plot description and save file names 
#' (if the results is saved to the outputs folder).
#' 
#' @param reshuffle_or_subset Either "reshuffle", or "subset" as a string. 
#' Instead of reshuffling a certain proportion, you can subset the metadata 
#' to not include that same proportion instead. This is an alternative form of
#' this analysis. Be aware that if this is set to "subset", anything that refers
#' to mismatch actually refers to the excluded cells, and anything that refers
#' to reshuffling refers to subsetting instead.
#' 
#' @param mismatch_props As a vector proportions between 0 and 1. To what degree
#' of mismatch should the cluster annotations be reshuffled to during the 
#' robusntess test. For example, c(0.1, 0.2, 0.3) would have the function 
#' compare the overlap between default CCI predictions and 10 % cluster 
#' annotation mismatch CCI predictions, as well as default and 20 % predictions,
#' and default and 30 % predictions.
#' 
#' Because a custom reshuffling function is used, the degree of mismatch and
#' degree of reshuffling are always the same.
#' 
#' At higher degrees of mismatch, some CCI inference methods produce errors 
#' because they can't find any significant predictions. As such it is 
#' recommended to sort your input vector from high-to-low proportions, so that
#' you can see if there are errors in the output early, and restart the run with
#' adjusted parameters. If you don't choose at least two proportions the
#' formatting of the code breaks.
#'
#' @param top_n As an integer. When ordering CCI predictions by relevance, where
#' should the cutoff for relevance be set? What "top_n" should be considered 
#' relevant? Almost every method has its own measurements for significance.
#' In order to compare methods, we must bridge the gaps between what each method
#' considers most relevant. To do so, the outputs of each method are ordered
#' according to that method's approach, and then the same top N ranks are 
#' for each method considered the "relevant". When we compare the overlap in
#' predicted CCIs, we compare use the overlaps of these top-ranks. 
#' 
#' @param number_seeds As an integer. Reshuffling has randomness to it, which
#' is why the reshuffling should be performed with multiple seeds and 
#' permutations to get a sample of the distribution of possible outcomes at a  
#' given mismatch proportion. As with any sample, the higher the "n" (in this 
#' case number_seeds), the more representative the sample is of the original 
#' distribution. 
#' 
#' @param methods_vector A vector of method names as strings. What methods 
#' should be included in the experiment? Make sure to use the name that LIANA++
#' uses for that method. Can only include methods in LIANA++.
#' 
#' @param liana_warnings Either TRUE, FALSE or "divert". Should the warnings 
#' from LIANA++, which are often repetitive and unhelpful, be either suppressed 
#' or alternatively diverted to the Logs folder? When these types of warning are
#' left in, they can often displace valuable warnings. Be careful with this 
#' setting, as suppressing warnings is obviously risky.
#' 
#' @param tag An additional tag that specifies the NATMI save files to this
#' iterator run.
#' 
#' @param save_results Either TRUE or FALSE. Should the plots and iterator
#' results objects be saved to the outputs folder with an automatically
#' generated file name? Recommended for keeping the results safe.
#' 
#' @param trial_run Is this a trial run of the iterator or serious results?
#' Takes a boolean. If this is a trial run, the save file names, logs and plot
#' captions will reflect this. 
#' 
#' @param cellchat_nperms As an integer. On some systems, CellChat is one of the
#' slower methods. It may be useful to set this parameter to 10 to speed up the 
#' analysis when running a test run. Unless you're using this function for a 
#' trial run, this should be left at its default, which is CellChat's own 
#' recommendation of 100.
#' 
#' @param outputs Which of the iterator results do you want in the final output?
#' By default, all the data generated is returned in a list, but in unusual 
#' scenarios you could subset this list with a vector passed to this argument. 
#' Each of the strings below is the name of a piece of data that can be returned
#' by adding it to the vector you pass here.
#' 
#' "top_ranks_overlap": A tibble of all the overlaps between original and
#' reshuffle-based predictions for each iteration, each reshuffling proportion
#' and each method.
#' 
#' "plot_box" and "plot_line": The two different plots that describe the 
#' top_ranks_overlap data.
#' 
#' "reshuffling_results": A bundled list of further information from the 
#' Iterator, including all the reshuffled cluster annotations, raw LIANA++
#' outputs, top-ranked CCIs for every condition, and so on.
#' 
#' "metadata": The metadata of the run. Save file names, parameters, run times 
#' and so on.
  
  
wrap_cluster_Iterator <-
  function(testdata,
           testdata_type,
           NATMI_tag,
           
           reshuffle_or_subset = "reshuffle",
           mismatch_props = c(seq(0.60, 0.05, -0.05)) ,
           top_n = 500,
           
           number_seeds = 10,
           methods_vector = c('call_connectome' ,
                              'call_natmi'      ,
                              'call_italk'      ,
                              'call_sca'        ,
                              'cellchat'        ,
                              'squidpy'),
           
           
           liana_warnings  = "divert",
           save_results    = TRUE,
           trial_run       = FALSE,
           
           
           
           cellchat_nperms = 100,
           
           outputs = c("top_ranks_overlap",
                       "plot_box",
                       "plot_line",
                       "reshuffling_results",
                       "metadata")) {
    
             
             
#------------------------------------------------------------------------------#
# 1.1 Generate Parameters  -----------------------------------------------------
{
  # In this segment, we take the user inputs in their simpler forms and turn 
  # into more functional parameters for our code
  
  
  
  # We generate number_ranks, a list with each methods name and the number of
  # top ranked CCIs to consider relevant for that method. Top_n is always the
  # same.
  number_ranks        <- as.list(rep(top_n, length(methods_vector)))
  names(number_ranks) <- methods_vector
  
  
  
  ## Formatting iterable structures
  {
    # The structures we create here can be iterated over in lapply's and maps.
    # They contain parameters we often need to iterate over and their formatting
    # is automatically carried over to the end result.
    
    
    # We format a named list of seeds, it contains as many seeds as the user
    # specified, from 1 to n, and each entry has an appropriate name, "Seed_n".
    seed_list <- as.list(1:number_seeds)
    
    names(seed_list) <-
      map(seed_list, function(seed) {
        # Name each element of seed_list appropriately
        str_glue("Seed_", seed)
        
      })
    
    

    # Convert mismatch_props to a list and name it.
    mismatch_props <- as.list(mismatch_props)
    
    names(mismatch_props) <- map(mismatch_props, function(prop) {
      # Name every Reshuffling proportion
      str_glue("Reshuffle_", as.character(prop * 100))
      
    })
    

    
    # Convert methods_vector to a list and name it
    methods_list        <- as.list(methods_vector)
    names(methods_list) <- methods_vector
    
  }
  

  
  # Not all testdata has consisten metadata and idents. We add a cluster_key 
  # column to the testdata, which is a set of cluster numerals counting from 
  # one upwards (cellchat throws an error if a cluster name is "0), and then 
  # make sure the Idents are equal to cluster_key.
  # Having the annotations be unifrom for every testdata makes our life a lot
  # easier in later steps.
  testdata@meta.data <- testdata@meta.data %>%
    mutate("cluster_key" = as.factor(as.numeric((Idents(testdata)))))
  
  Idents(testdata) <-  testdata@meta.data$cluster_key
  
  # If the Idents of the testdata are not equal any of the metadata columns, we 
  # throw an error, because that means the previous step didn't work.
  if(is.null(liana:::.get_ident(testdata))) {
    stop(str_glue("There is no column in the metadata of the seurat object ",
                  "that is equal to the seurat object's idents"))
    
  }
  
  
  
  # Format the Sys.time() of the run to not contain characters that are bad to
  # have in save file names. We will later use this to tag file names and
  # plots so they can be grouped according to run, and all have unique names.
  time_of_run <-  Sys.time() %>%
    as.character()    %>%
    gsub(':', '-', .) %>% # save files can't have colons
    gsub(' ', '_', .) %>% # save files shouldn't have spaces
    str_sub(1 , nchar(.) - 3) # the code never runs in under a minute, so the
                              # number of seconds isn't valuable information.
  
  
  # If necessary we generate a filepath to save LIANA++ logs under.
  if (liana_warnings == "divert") {
    warning_logfile <-
      clust_auto_file_Name(
        prefix = "Outputs/Cluster_Reshuffling/Logs/",
        suffix = ".txt",
        
        reshuffle_or_subset = reshuffle_or_subset,
        testdata_type       = testdata_type,
        number_ranks        = number_ranks,
        time_of_run         = time_of_run,
        trial_run           = trial_run)
    
  }
  
  # Generate the filepaths to save the data under, if necessary.
  # CR stands for Cluster Reshuffling.
  if (save_results == TRUE) {
    
    box_plot_png_name <-
      clust_auto_file_Name(
        prefix = "Boxplot_CR_",
        suffix = ".png",
        
        reshuffle_or_subset = reshuffle_or_subset,
        testdata_type       = testdata_type,
        number_ranks        = number_ranks,
        time_of_run         = time_of_run,
        trial_run           = trial_run)
    
    line_plot_png_name <-
      clust_auto_file_Name(
        prefix = "Lineplot_CR_",
        suffix = ".png",
        
        reshuffle_or_subset = reshuffle_or_subset,
        testdata_type       = testdata_type,
        number_ranks        = number_ranks,
        time_of_run         = time_of_run,
        trial_run           = trial_run)
    
    iterator_results_save_path <- 
      clust_auto_file_Name(
        prefix = "Outputs/Cluster_Reshuffling/Iterator_Results_CR_",
        suffix = ".RData",
        
        reshuffle_or_subset = reshuffle_or_subset,
        testdata_type       = testdata_type,
        number_ranks        = number_ranks,
        time_of_run         = time_of_run,
        trial_run           = trial_run)
    
  }
  
  

} 


    
#------------------------------------------------------------------------------#
# 1.2 Create Reshuffled Cluster Annotations ------------------------------------
{
  # We reach for the shuffler or subsetter function depending on user input.
  if (reshuffle_or_subset == "reshuffle") {
    
    # Let the user follow along in the console
    print_Title("1. Creating Reshuffled Cluster Annotations",
                space_after = 0, 
                super = TRUE)
    
    
    
    # For every mismatch proportion and for every seed, generate a meta.data df 
    # with reshuffled cluster annotations.
    reshuffled_clusters <- lapply(
      mismatch_props,
      wrap_Shuffler,
      seed_list = seed_list,
      metadata  = testdata@meta.data) %>%
      # Add the default metadata to the reshuffled cluster for completeness
      append(list("Reshuffle_0" = testdata@meta.data), .)
    
    
    
  } else if (reshuffle_or_subset == "subset") {
    
    # Let the user follow along in the console
    print_Title("1. Subsetting Cell Clusters",
                space_after = 0, 
                super = TRUE)
    
    
    # define the proportion of the cluster to remove 
    removal_props <- mismatch_props
    
    # For every removal proportion and for every seed, generate a meta.data df 
    # with subset cluster annotations.
    reshuffled_clusters <- lapply(
      removal_props,   # mismatch props is actually a subset props
      wrap_Subsetter,
      seed_list = seed_list,
      metadata  = testdata@meta.data) %>%
      # Add the default metadata to the reshuffled cluster for completeness
      append(list("Reshuffle_0" = testdata@meta.data), .)
    
    
    
  }
   
  

  
  
  
}


    
#------------------------------------------------------------------------------#
# 1.3 Run LIANA ----------------------------------------------------------------
{
  # iterate_liana_wrap() runs LIANA++ on testdata with every cluster annotation
  # we created (so again, default, then for every seed and every mismatch prop)
  liana_results <-
    iterate_liana_wrap(
      seed_list           = seed_list,
      mismatch_props      = mismatch_props,
      reshuffled_clusters = reshuffled_clusters,
      testdata            = testdata,
      methods_vector      = methods_vector,
      
      liana_warnings   = liana_warnings,
      warning_logfile  = warning_logfile,
      
      tag = NATMI_tag,
      
      # Add parameters for liana_wrap that are passed as (...)
      expr_prop = 0.1,
      cellchat.params   = list(nboot     = cellchat_nperms, 
                               expr_prop = 0.1,
                               thresh    = 1)
    )  
  
  
  # The output has the runtime tacked on to the end of it, we separate that here
  runtime <- liana_results$runtime
  liana_results$runtime <- NULL
  
  # The easiest structure for iteration is not the best for accessing the data
  # We restructure the hierarchy of the results through multiple transpositions
  liana_results <- liana_results %>%
    map_depth(., .depth = 1, transpose) %>%
    map_depth(., .depth = 0, transpose) %>%
    
    # Correct Cellchat Errors where nothing is significant
    map_depth(., .depth = 3, function(result) {
      
      if (is_tibble(result) == FALSE) {
        if (result$message == str_glue("No significant signaling interactions ",
                                       "are inferred based on the input!")) {
          
          # replace the error with an empty tibble, non significance means
          # 0 predicted interactions, not an error.
          return(liana_results$Reshuffle_0$Seed_1$cellchat[0, ])
          
        }
        
      }
      
      return(result)
      
    }) %>%
    
    # If an error occurs in LIANA++ it returns it instead of an output tibble
    # Here we check if any errors (= non-tibbles) were returned.
    map_depth(., .depth = 3, function(result) {
      
      # If LIANA++ didn't return a tibble here, save the results so far and 
      # metadata. The rest of the code won't work from here on out.
      # LIANA is the most likely part of the code for errors to be generated, 
      # which is why we check for it here and salvage what's left if something
      # went wrong.
      if(is_tibble(result) == FALSE) {
        
        # summarize the metadata
        metadata <- clust_summarise_Metadata(
          seed_list        = seed_list,
          mismatch_props   = mismatch_props,
          methods_list     = methods_list,
          
          testdata_type       = testdata_type,
          reshuffle_or_subset = reshuffle_or_subset,
          number_ranks        = number_ranks,
          
          cellchat_nperms  = cellchat_nperms,
          outputs          = outputs,
          
          liana_warnings   = liana_warnings,
          save_results     = save_results,
          trial_run        = trial_run,
          
          runtime     = runtime,
          time_of_run = time_of_run,
          
          warning_logfile    = warning_logfile,
          line_plot_png_name = line_plot_png_name,
          box_plot_png_name  = box_plot_png_name,
          iterator_results_save_path = iterator_results_save_path
        )
        
        # Mash all the results so far together
        error_results <- list("metadata" = metadata,
                              "liana_results" = liana_results,
                              "reshuffled_clusters" = reshuffled_clusters)
        
        # Save what we can, with an error warning in the file name
        # We subset the original file name to paste in the ERROR part
        save(error_results, 
             file = str_glue(str_sub(iterator_results_save_path,
                                     1, 
                                     nchar(iterator_results_save_path) - 6),
                             "_ERROR.RData"))
        
        # Explain the situation to the user
        stop(str_glue("An error occured in one of the LIANA methods. ",
                      "Instead of an output tibble, LIANA returned: \n",
                      as.character(result), 
                      "\n\n The mismatch proportion may be too high to return ", 
                      "any significant LR Interactions, this makes some ",
                      "methods, such as CellChat crash."))
        
      } else {
        
        # if nothing is wrong do nothing and proceed.
        return(result)
        
      }
      
    })
  
  
}



#------------------------------------------------------------------------------#
# 1.4 Compare Top-Ranked Predictions -------------------------------------------
{
  # In this segment we determine the top-ranked-interactions for every condition
  # and compare them to the baseline predictions.
  
  
  # Segment header so the user can follow along in con
  print_Title(str_glue("4. Calculate Overlap Between Default and Reshuffled ",
                       "CCI predictions"),
              super = TRUE)
  
  
  
  top_ranks <-
    # We extract the top ranked interactions of all our LIANA++ results
    map(methods_list, function(method) {
      # clust_get_top_ranks needs to know what method its working on, so we
      # map by method
      top_ranks_for_method <- liana_results[[method]] %>%
        # Now we use map_depth to map down to the tibbles in the nested list
        # and apply our function to them
        map_depth(
          .,
          .depth = 2,
          clust_get_top_ranks,
          method = method,
          top_n = number_ranks[[method]],
          with_ties = TRUE
        )
      
    })  %>%
    # Now that we have the top_ranks in a nested list we map to the lowest level
    # and add an LR_Pair and LR_ID column to the tibbles, which we need to 
    # compare overlaps.
    map_depth(., .depth = 3, format_top_ranks)
  
  
  
  
  # We map over all the top-ranks-filled tibbles and compare the overlap
  # between the original predictions and various later ones using the 
  # rank_overlap utils function.
  top_ranks_overlap <- map(methods_list, function(method) {
    
    # We map by method because each methods base line predictions are different
    overlaps_for_method <-  top_ranks[[method]] %>%
      map_depth(.,
                .depth = 2,
                rank_overlap,
                main_ranks = top_ranks[[method]]$Reshuffle_0[[1]],
                verbose = FALSE,
                expect_same_size = FALSE)
    
  }) %>%
  # Turning top_ranks_overlap into one concise tibble (also helps plotting)
  # This code won't work if only one mismatch_prop was selected.
    as_tibble()                             %>%
    mutate(Mismatch = c(0, mismatch_props)) %>%
    unnest(cols = all_of(methods_vector))   %>%
    unnest(cols = everything())             %>%
    relocate(Mismatch)                      %>%
    pivot_longer(cols = !(starts_with("Mismatch")), 
                 names_to = "Method")       %>%
    arrange(Mismatch)                       %>%
    arrange(Method)                         %>%
    rename("Overlap" = value)
  
  
  # print a sample of the output for the user
  print(head(top_ranks_overlap, 15))
  
  
}


    
#------------------------------------------------------------------------------#
# 1.5 Plotting of Results ------------------------------------------------------
{
  
  # Here we visualize the overlap between top ranked CCI predictions as they
  # change with the reshuffling of cluster annotations. 
  # once as a boxplot, and once as a scatter/line plot.
  
  
  # We reformat the collated_top_ranks_overlap tibble so its more suitable for
  # plotting
  
  # The Reshuffling proportion and overlap are clearer in percentage
  # Rename the methods from the LIANA++ internal string to their official name
  tr_overlap_for_plot <- top_ranks_overlap %>%
    as.data.frame()                        %>%
    mutate(Mismatch = Mismatch * 100)      %>% # proportion to percent
    mutate(Overlap  = Overlap  * 100)      %>% # proportion to percent
    as_tibble()                            %>%
    mutate("Method" = recode(Method,
                             "call_connectome" = "Connectome",
                             "squidpy"         = "CellPhoneDB",
                             "call_natmi"      = "NATMI",
                             "call_italk"      = "LogFC Product",
                             "call_sca"        = "SingleCellSignalR",
                             "cellchat"        = "CellChat")) # renaming
  
  
  
  
  # To directly be able to associate the box plot with the settings that
  # produced it, we automatically generate a plot description
  plotting_caption <-
    clust_plot_Description(
      mismatch_props      = mismatch_props,
      trial_run           = trial_run,
      testdata_type       = testdata_type,
      reshuffle_or_subset = reshuffle_or_subset, 
      seed_list           = seed_list,
      number_ranks        = number_ranks,
      time_of_run         = time_of_run
    )
  
  
  if (reshuffle_or_subset == "reshuffle") {
    
    plot_line <- 
      overlap_line_Plot(tr_overlap_for_plot,
                        plotting_caption,
                        x_axis_var  = "Mismatch",
                        x_axis_name = "Mismatch of Cluster Annotations [%]")
    
    plot_box <- 
      overlap_box_Plot(tr_overlap_for_plot,
                       plotting_caption,
                       x_axis_var  = "Mismatch",
                       x_axis_name = "Mismatch of Cluster Annotations [%]")
    
    
  } else if (reshuffle_or_subset == "subset") {
    
    plot_line <- 
      overlap_line_Plot(tr_overlap_for_plot,
                        plotting_caption,
                        x_axis_var  = "Mismatch",
                        x_axis_name = "% of Cells Removed per Cluster")
    
    plot_box <- 
      overlap_box_Plot(tr_overlap_for_plot,
                       plotting_caption,
                       x_axis_var  = "Mismatch",
                       x_axis_name = "% of Cells Removed per Cluster")
    
    
  }
  
  # Generate our plots with functions.
  
  # Print out visualizations
  print(plot_line)
  print(plot_box)
  
  
  
  
  # Removing Clutter
  rm(tr_overlap_for_plot, plotting_caption)
  
  
}



#------------------------------------------------------------------------------#
# 1.6 Capturing Script Metadata ------------------------------------------------
{
  # In this segment, we summarize the metadata of this Iterator run. When 
  # troubleshooting or reproducing results, this information will be useful to 
  # the user.
  
  
  # Segment header for the user to follow along in the console
  if (save_results == TRUE) {
    print_Title(str_glue("5. Summarizing Metadata and Saving Results"),
                super = TRUE)
    
  } else {
    print_Title(str_glue("5. Summarizing Metadata"),
                super = TRUE)
    
  }
  

  
  # We already have the points in time from the LIANA iterator that refer to 
  # when the most time consuming steps took place. We use this named list pf
  # time points to create a succinct and informative runtime overview.
  runtime <- calculate_Runtime(runtime)
  
  
  
  # We then summarise the above information and more metadata into a single
  # object
  metadata <- clust_summarise_Metadata(
    seed_list        = seed_list,
    mismatch_props   = mismatch_props,
    methods_list     = methods_list,
    
    testdata_type       = testdata_type,
    reshuffle_or_subset = reshuffle_or_subset,
    number_ranks        = number_ranks,
    
    cellchat_nperms  = cellchat_nperms,
    outputs          = outputs,
    
    liana_warnings   = liana_warnings,
    save_results     = save_results,
    trial_run        = trial_run,
    
    runtime     = runtime,
    time_of_run = time_of_run,
    
    warning_logfile    = warning_logfile,
    line_plot_png_name = line_plot_png_name,
    box_plot_png_name  = box_plot_png_name,
    iterator_results_save_path = iterator_results_save_path
  )
  
  
  
  # Now that these objects are stored in the metadata object, we can remove
  # this clutter from the environment.
  rm(runtime, seed_list, mismatch_props, number_ranks, methods_list,
     cellchat_nperms, methods_vector, number_seeds, testdata_type, 
     time_of_run, trial_run)
  
  
}



#------------------------------------------------------------------------------#
# 1.7 Packaging Results to return them -----------------------------------------
{
  # In this segment we bundle up all the findings into one returnable object
  
  
  # Bundle up the non-essential cluster annotations, liana results and so on.
  reshuffling_results <- list("testdata"      = testdata,
                              "reshuffled_clusters" = reshuffled_clusters,
                              "liana_results" = liana_results,
                              "top_ranks"     = top_ranks)
  
  # Remove associated clutter
  rm(reshuffled_clusters,testdata,liana_results, top_ranks)
  
  
  
  # Bundle up all our findings into one object
  iterator_results <-
    list(
      "top_ranks_overlap"  = top_ranks_overlap,
      "plot_box"  = plot_box,
      "plot_line" = plot_line,
      "reshuffling_results" = reshuffling_results,
      "metadata"  = metadata
    )
  
  
  
  # Filter our results by the outputs the user wants to retrieve
  # Usually´all the data is requested so this step doesn't change anything.
  iterator_results <- iterator_results[outputs]
  
  # Remove associated clutter
  rm(top_ranks_overlap, reshuffling_results, metadata)
  
  
}



#------------------------------------------------------------------------------#
# 1.8 Saving Results -----------------------------------------------------------
{
  # In this segment we save the results and let the user know where to find all
  # the files that were generated.
  
  
  # clust_save_Results() prints where files were saved to, we expand on that
  # here by letting the user know where the LIANA warning file is stored.
  if(liana_warnings == "divert") {
    
    # let the user know where to find the log
    cat(str_wrap(str_glue("LIANA warnings saved at ~/", 
                          warning_logfile, "."), width = 60), " \n\n")
    
    
  }
  
  
  
  
  # In this segment we save the plots and environment to the outputs folder,
  # if it's specified in by the user. This also prints the respective file paths
  # to the console.
  if (save_results == TRUE) {
    clust_save_Results(
      plot_box  = plot_box,
      plot_line = plot_line,
      iterator_results = iterator_results,
      
      line_plot_png_name = line_plot_png_name,
      box_plot_png_name  = box_plot_png_name,
      iterator_results_save_path = iterator_results_save_path)
    
  }
  
  
}
  

    
# And now return our results to the user
return(iterator_results)

      
} # end of function

