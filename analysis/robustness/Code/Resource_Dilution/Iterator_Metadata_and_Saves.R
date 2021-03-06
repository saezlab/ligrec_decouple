#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # This is the Iterator_Metadata_and_Saves.R script.

  # The goal of this script is do define the more meta level functions in the
  # RD_Robustness_Iterator.R script that don't directly process data, but
  # instead are relevant for saving files and summarizing metadata.

}



#------------------------------------------------------------------------------#
# 1. Defining Functions---------------------------------------------------------


# summarise_Metadata()
{
  #' Summarizes the metadate relevant for the Robustness Iterator
  #'
  #' @description This function summarizes all the iterator parameters, file
  #' names that were used, runtime data, and more into one metadata object.
  #'
  #'
  #' @return Returns a compiled list of metadata, parameters and save file
  #' locations (if files were saved to the computer).



  summarise_Metadata <- function(number_seeds,
                                 testdata_type,
                                 feature_type,
                                 preserve_topology,
                                 modify_baseline,
                                 dilution_props,
                                 top_n,
                                 methods_vector,

                                 liana_warnings,

                                 cellchat_nperms,
                                 bundled_outputs,
                                 master_outputs,

                                 save_results,
                                 trial_run,

                                 runtime,
                                 time_of_run,

                                 warning_logfile,
                                 box_plot_rds_name,
                                 line_plot_rds_name,
                                 iterator_results_save_path) {

    # Summarize the metadata parameters
    meta_params <- list(
      "time_of_run"  = time_of_run,
      "save_results" = save_results,
      "trial_run"    = trial_run
    )

    # summarise all the parameters from wrap_resource_Iterator()
    dilution_params <- list(
      "number_seeds"      = number_seeds,
      "testdata_type"     = testdata_type,
      "feature_type"      = feature_type,
      "modify_baseline"   = modify_baseline,
      "preserve_topology" = preserve_topology,
      "dilution_props"    = dilution_props,
      "top_n"             = top_n,
      "methods_vector"    = methods_vector,

      "liana_warnings"    = liana_warnings,

      "cellchat_nperms"   = cellchat_nperms,
      "bundled_outputs"   = bundled_outputs,
      "master_outputs"    = master_outputs
    )

    # Put all the parameters in a list
    metadata <- list(
      "runtime"         = runtime,
      "dilution_params" = dilution_params,
      "meta_params"     = meta_params,
      "sessionInfo"     = sessionInfo()
    )


    # If the results were saved, tack the file names that the saves are under
    # onto the end of metadata.
    if (save_results == TRUE) {
      # Generate the filepaths data was saved under.
      # RD stands for Resource Dilution.
      metadata[["box_plot_rds_name"]] <-
        box_plot_rds_name

      metadata[["line_plot_rds_name"]] <-
        line_plot_rds_name

      metadata[["iterator_results_save_path"]] <-
        iterator_results_save_path
    }



    # If a warnings log was created for resource_Robustness(), append the file
    # name of the log to the metadata.
    if (liana_warnings == "divert") {

      metadata[["warning_logfile"]] <-
        warning_logfile

    }


    # return the metadata.
    return(metadata)

  } # end of function



}


# save_Results()
{
  #' Saves its three arguments to custom filepaths
  #'
  #' @param plot_box Takes the boxplot generated by overlap_box_Plot as an
  #' input. Saves it to the outputs folder under a descriptive name.
  #'
  #' @param plot_line Takes the lineplot generated by overlap_line_Plot as an
  #' input. Saves it to the outputs folder under a descriptive name.
  #'
  #' @param iterator_results Takes the list of results from the iterator, saves
  #' them to a descriptive file name in the outputs folder.
  #'
  #' @param box_plot_rds_name Name for the boxplot, ending in ".rds".
  #' @param line_plot_rds_name Name for the lineplot, ending in ".rds".
  #' @param iterator_results_save_path Filepath to save the iterator results
  #' under.

  save_Results <- function(plot_box,
                           plot_line,
                           iterator_results,

                           box_plot_rds_name,
                           line_plot_rds_name,
                           iterator_results_save_path) {


    # Save both plots
    saveRDS(
       object = plot_box,
       file = str_glue(
           "analysis/robustness/Outputs/Resource_Dilution/",
           box_plot_rds_name
       )
    )

    saveRDS(
        object = plot_line,
        file = str_glue(
            "analysis/robustness/Outputs/Resource_Dilution/",
            line_plot_rds_name
        )
    )

    # Save R environment and all the results within it
    save(iterator_results, file = iterator_results_save_path)




    # Let the user know where everything was stored.
    cat(str_wrap(str_glue("Box Plot saved at ~/analysis/robustness/Outputs/Resource_Dilution/",
                   box_plot_rds_name, "."), width = 60), "\n\n")

    cat(str_wrap(str_glue("Line Plot saved at ~/analysis/robustness/Outputs/Resource_Dilution/",
                   line_plot_rds_name, "."), width = 60), "\n\n")

    cat(str_wrap(str_glue("Iterator Results saved at ~/",
                   iterator_results_save_path, "."), width = 60), "\n\n")

  }  # end of function
}


# auto_file_Name()
{
  #' Automatically generates a file name or file path
  #'
  #' @param prefix As a char. What should the file name start with? It could be
  #' a folder to make it a file path, such as "Outputs/", or any other tag, or
  #' "".
  #'
  #' @param suffix As a char. What should the file name end with? It should be a
  #' file extension such as ".txt" or ".RData" at minimum, but it could also be
  #' more, such as "report.txt".
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
  #' @return A file name that starts with the prefix, ends with the suffix and
  #' contains a bunch of parameter tags in between. This way the user can
  #' identify the save file by the parameters it was set up with.


  auto_file_Name <- function(prefix,
                             suffix,

                             trial_run,
                             preserve_topology,
                             modify_baseline,
                             testdata_type,
                             feature_type,
                             number_ranks,
                             time_of_run) {

    # We define individual comments related to relevant parameters and then
    # string them all together for the save file name.

    # If this is a trial run, mark the save files as such
    if (trial_run == FALSE) {

      test_run_comment <- ""

    } else if (trial_run == TRUE) {

      test_run_comment <- "TRIAL_RUN_"

    }

    # If the baseline could be modified, mark the save files as such
    if (modify_baseline == FALSE) {

      baseline_comment <- ""

    } else if (modify_baseline == TRUE) {

      baseline_comment <- "mod_baseline_"

    }


    # How was the topology of diluted resources handled? Make an appropiate
    # comment.
    if (preserve_topology == FALSE) {
      topology_comment <- "rand_topo_"

    } else if (preserve_topology == TRUE) {
      topology_comment <- "pres_topo_"

    }


    # ´What testdata_type was extracted and used with resource_Robustness?
    testdata_comment <-
      str_glue(testdata_type, "_")


    # What feature_type was dilution specified to occur with?
    feature_type_comment <-
      str_glue(feature_type, "_")


    # Make a comment out of the median top_n that was considered top_ranked.
    # As a note, it could be a different number per method, but often it will
    # be the same number for every method.
    top_ranks_comment <-
      str_glue("top", median(unlist(number_ranks)), "_", )


    # Mash all the comments together with the suffix and prefix to create our
    # custom and hopefully informative file names
    auto_file_name <-
      str_glue(
        prefix,
        test_run_comment,
        testdata_comment,
        baseline_comment,
        topology_comment,
        feature_type_comment,
        top_ranks_comment,
        time_of_run,
        suffix
      )

    return(auto_file_name)
  }  # end of function

}

