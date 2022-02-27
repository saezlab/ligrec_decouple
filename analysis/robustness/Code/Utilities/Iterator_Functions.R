#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # The Iterator_Functions.R script. As a utilities script, it defines functions
  # that are used in both resource dilution and cluster reshuffling. The
  # functions here are grouped because they all play a role in running the
  # iterator for each of these respective processes.

}



#------------------------------------------------------------------------------#
# 1. Define Functions ----------------------------------------------------------

# get_top_ranks()
{
  #' Get the top n ranked items of a method from the tibble liana wrapper or
  #' call_x results
  #'
  #' @param data_set The tibble output by the liana wrapper function or call_x
  #' functions.
  #' @param method The method for which you would like to extract the top ranked
  #' interactions, as a string.
  #' @param top_n The number of items to return, returns items ranked #1 to #n.
  #'
  #' @return Returns the tibble input cut down to the top n highest ranked
  #' interactions.

  get_top_ranks <-
    function(data_set, method, top_n) {
      # generate a list that describes how to rank in each method to get what the
      # method considers best
      rank_spec_list <-
        list(
          "cellchat" = list(method_score = "prob",
                            descending_order =  TRUE),

          "call_connectome" = list(method_score = "weight_sc",
                                   descending_order =  TRUE),

          "call_italk" = list(method_score = "logfc_comb",
                              descending_order =  TRUE),

          "call_natmi" = list(method_score = "edge_specificity",
                              descending_order =  TRUE),

          "call_sca" = list(method_score = "LRscore",
                            descending_order =  TRUE),

          "cellphonedb" = list(method_score = "lr.mean",
                               descending_order =  TRUE),

          "cytotalk" = list(method_score = "crosstalk_score",
                            descending_order = TRUE),

          "logfc" = list(method_score = "logfc_comb",
                         descending_order =  TRUE)
        )



      # If its a p-value method, the code will go into this if statement
      if (rank_spec_list[[method]]$descending_order == FALSE) {
        topranks <-
          slice_min(
            data_set,
            n = top_n,
            order_by = !!sym(rank_spec_list[[method]]$method_score),
            with_ties = TRUE
          )

        # order_by is assigned the criterion dictated by rank_spec_list


      } else {
        # if it's not one of the p-value methods

        topranks <-
          slice_max(
            data_set,
            n = top_n,
            order_by = !!sym(rank_spec_list[[method]]$method_score),
            with_ties = TRUE
          )

        # order_by is assigned the criterion dictated by rank_spec_list


      }

      return(topranks)

    } #end of function
}

# extract_Testdata()
{
  #' Helper function that gets a specific seurat object from the outputs folder
  #'
  #' @param testdata_type As a string. Which testdata should be retrieved?
  #' Either "seurat_pbmc" or "liana_test". Seurat_pbmc is the data set used in
  #' the seurat tutorial, while liana_test is the testdata that comes with
  #' LIANA++, and is a small subset of seurat_pbmc.
  #'
  #' @return A seurat object loaded from the outputs folder or liana package.


  extract_Testdata <- function(testdata_type) {

    # Get seurat or liana test data
    if (testdata_type == "seurat_pbmc") {

      # Read testdata from outputs
      testdata <- readRDS(file = "analysis/robustness/Data/pbmc3k_final.rds")

    } else if (testdata_type == "liana_test") {

      # Where is the liana testdata located?
      liana_path <- system.file(package = 'liana')
      # Read the testdata from its location.
      testdata <-
        readRDS(file.path(liana_path, "testdata", "input", "testdata.rds"))

      # removing superfluous values
      rm(liana_path)


    } else {

      # error if its not one of the supported data sets
      stop("Testdata name not recognized!")

    }


    # Return the seurat.
    return(testdata)

  } # end of function
}


#' Takes get_n_top_ranks outputs that have an LR_ID and determines their
#' Jaccard Index
#'
#' @description An LR_ID is a unique identifier for a CC Interaction. It is
#' simply the concatenated source name, target name, ligand name and receptor
#' name. A CCI is fully characterized by its LR_ID, which is why they can be
#' used when comparing two tibbles with top ranked interactions.
#'
#' @param main_ranks A tibble of of top ranked interactions (Usually no dilution,
#'  i.e  ground truth interactions)
#' @param comparison_ranks A tibble of top ranked interactions
#' @param verbose Should the function describe the JI to you or not?
#'
#' @return The JI (0-1) between the two input frames in contents of the
#' LR_ID column, as well as an optional print statement that gives more
#' detail.
rank_overlap <- function(main_ranks,
                         comparison_ranks,
                         verbose = TRUE) {

  # calculate overlap between LR_IDs
  LRID_intersect <- sum(comparison_ranks$LR_ID %in% main_ranks$LR_ID)
  LRID_union     <- c(comparison_ranks$LR_ID, main_ranks$LR_ID) %>%
    unique() %>%
    length()

  jaccard_index <- LRID_intersect / LRID_union


  # describe the output to the user
  if (verbose == verbose) {
    print(str_glue(""))
    cat(
      str_wrap(
        str_glue(
          "The main ranking and the comparison ranking have ",
          as.character(LRID_intersect),
          " LR_IDs in common. This corresponds to a jaccard index of ",
          round(jaccard_index, 2),
          "."),
        width = 60), "\n")
  }

  return(jaccard_index)

}

#' Takes get_n_top_ranks outputs that have an LR_ID and determines their
#' TPR
#'
#' @description An LR_ID is a unique identifier for a CC Interaction. It is
#' simply the concatenated source name, target name, ligand name and receptor
#' name. A CCI is fully characterized by its LR_ID, which is why they can be
#' used when comparing two tibbles with top ranked interactions.
#'
#' @param main_ranks A tibble of of top ranked interactions (Usually no dilution,
#'  i.e  ground truth interactions)
#' @param comparison_ranks A tibble of top ranked interactions
#' @param verbose Should the function describe the calc or not
tpr_overlap <- function(main_ranks, # p_ranks (only condition positives)
                        comparison_ranks, # all test ranks
                        verbose = TRUE) {

  # calculate TPR
  tp <- sum(comparison_ranks$LR_ID %in% main_ranks$LR_ID)
  cp <- length(comparison_ranks$LR_ID)
  tpr <- tp / cp # cp = all condition positives

  return(tpr)

}




# calculate_Runtime()
{
  #' Converts a list of checkpoint names and associated system times into a
  #' convenient tibble that highlights the time that passed between the
  #' checkpoints
  #'
  #' @description Takes the runtime output of resource_Robustness(), which is a
  #' named list of checkpoints in time, and creates a  output tibble that has
  #' the time between each checkpoint in it and the time elapsed up until that
  #' checkpoint.
  #'
  #' @param runtime A list of checkpoints from the Iterator and the times at
  #' which they were reached.
  #'
  #' @return A tibble giving an overview of the runtime of a piece of code,
  #' typically of the Iterator..


  calculate_Runtime <- function(runtime) {

    # save the names of the time-points for later
    runtime_labels <- names(runtime)

    # convert run time to numeric so we can perform arithmetic operations on
    # them. In this case we need it for subtractions, to calculate the duration
    # between checkpoints
    runtime_numeric <- as.numeric(runtime)

    # We calculate the passage of time between checkpoints in the
    # resource_Robustness().
    # Step duration is the duration of a step between neighboring checkpoints.
    # Time elapsed is the duration between the completion of a step and the
    # start of the script.

    step_duration <- c(0) # No time has passed when the script is initialized.
    time_elapsed  <- c(0) # No time has passed when the script is initialized.

    # starting with the second index of runtime_numeric until the last index
    for (i in 2:length(runtime_numeric)) {

      # subtract the preceding checkpoint from the checkpoint at i, this is the
      # amount of time that passed between these two checkpoints
      step_duration <- c(step_duration,
                         runtime_numeric[[i]] - runtime_numeric[[i-1]])

      # subtract the very first checkpoint from the checkpoint at i, this is all
      # the time that has elapsed up until now.
      time_elapsed  <- c(time_elapsed,
                         runtime_numeric[[i]] - runtime_numeric[[1]])

    }

    # Turn seconds into time periods using lubridate and round for simplicity
    # Time periods are HH:MM:SS, which is earier to understand than just values
    # in seconds.
    step_duration <- round(seconds_to_period(step_duration))
    time_elapsed  <- round(seconds_to_period(time_elapsed))


    # summarize all the runtime data in a tibble
    runtime <- runtime %>%
      as_tibble_col() %>%
      unnest(cols = c(value)) %>%
      rename("Start Time" = "value") %>%
      add_column("Step Name" = runtime_labels, .before = 1) %>%
      add_column("Step Duration" = step_duration) %>%
      add_column("Time Elapsed" = time_elapsed)


    # Get rid of clutter in the environment
    rm(runtime_numeric,
       step_duration,
       time_elapsed,
       runtime_labels,
       i)

    # return the runtime tibble
    return(runtime)
  }
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


#' Function to calculate FPR from robustness results
#'
#' @param result_path path to the .rdata with the results from all runs.
#'
#' @return a tibble method, shuffle, and fpr
calculate_fpr <- function(result_path, analysis_focus){
  # load data
  load(result_path)

  # extract top ranks from all methods
  if(analysis_focus == "cluster"){
    top_ranks <- iterator_results$reshuffling_results$top_ranks
    test_ranks <- iterator_results$reshuffling_results$liana_results
  } else if(analysis_focus == "resource"){
    top_ranks <- iterator_results$collated_robustness_results$top_ranks_OP
    test_ranks <- iterator_results$collated_robustness_results$liana_results_OP
  }

  # Count N (i.e. number of Condition negatives)
  cn_list <- imap(top_ranks, function(top_res, method_name){

    # establish ground truth
    if(analysis_focus == "cluster"){
      all_gt_ranks <- test_ranks[[method_name]]$Reshuffle_0$Seed_1
      top_gt_ranks <- top_res$Reshuffle_0$Seed_1

    } else if(analysis_focus == "resource"){
      all_gt_ranks <- test_ranks[[method_name]]$OmniPath_0_Seed_1
      top_gt_ranks <- top_res$OmniPath_0_Seed_1
    }

    nrow(anti_join(all_gt_ranks,
                   top_gt_ranks,
                   by = c("source", "target",
                          "ligand", "receptor"))
    )
  })


  # calculate number of FPs in regards to 0 reshuffling
  # (We only need the top ranks from both Condition and Test)
  # We map test top_ranks
  res_to_fpr <- imap(top_ranks, function(method_res, method_name){

    # establish ground truth
    if(analysis_focus == "cluster"){
      top_gt_res <- top_ranks[[method_name]]$Reshuffle_0$Seed_1

    } else if(analysis_focus == "resource"){
      top_gt_res <- top_ranks[[method_name]]$OmniPath_0_Seed_1
    }

    method_res %>%
      enframe(name='shuffle', value = 'results') %>%
      {`if`(analysis_focus == 'cluster', unnest(., results), .)} %>%
      mutate(n_fp = map_dbl(results, function(res){
        count_fp(top_test_ranks = res,
                 top_gt_ranks = top_gt_res)
      })
      ) %>%
      mutate(fpr = n_fp / cn_list[[method_name]]) %>%
      select(-results)
  }) %>%
    enframe(name = 'Method') %>%
    unnest(value) %>%
    # Remove seed from resource-focused results
    mutate(shuffle = gsub("_Seed_.", "", shuffle)) %>%
    separate(shuffle, into = c("x", "value"), remove = FALSE) %>%
    select(-x) %>%
    mutate(value = as.numeric(value)) %>%
    mutate(Method = recode_methods(Method))

  return(res_to_fpr)
}

#' Helper function to count the number of FPs
count_fp <- function(top_test_ranks,
                     top_gt_ranks){
  n_fp <- sum(!(top_test_ranks$LR_ID %in% top_gt_ranks$LR_ID))
  return(n_fp)
}

