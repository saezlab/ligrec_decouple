#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # The User_Outputs_and_Plots.R script. As a utilities script, it defines
  # functions that are used in both resource dilution and cluster reshuffling.
  # The functions here are grouped because they all serve to create outputs for
  # the user about the Iterator. For example, print_Title is a helper function
  # to display in the console what stage the code is running at, while the
  # plotting functions give the user insight into the end results of the
  # Iterator.

}



#------------------------------------------------------------------------------#
# 1. Define Functions for Plotting Outputs -------------------------------------

# overlap_line_Plot()
{

  #' Draw a line/scatter plot from the collated top_ranks_overlap data
  #'
  #' @param x_axis_var What variable should the x-axis of the plot be? In our
  #' two setups, this is usually the percentage of dilution for a resource or
  #' the amount of mismatch in the reshuffled cluster annotations.
  #'
  #' @param x_axis_name The x-axis label that will be used in the plot.
  #'
  #' @param top_ranks_overlap As a tibble, formatted with a column for the
  #' x-axis, a column for the method name, and a column for the overlap of top
  #' ranked CCIs for that method for that x-axis value.
  #'
  #' @param plotting_caption A caption for your plot, as a string. Potentially
  #' could be an auto_plot_Description() output.
  #'
  #' @return A ggplot of the top_ranks_overlap as a line and scatter plot,
  #' complete with a caption.


  overlap_line_Plot <- function(top_ranks_overlap,
                                plotting_caption,
                                x_axis_var,
                                x_axis_name) {

    plot_line <-
      ggplot(data = top_ranks_overlap, aes(.data[[x_axis_var]],
                                           Overlap,
                                           group = Method,  # for colors
                                           color = Method)) +

      # add transparency because some points may overlap
      geom_point(alpha = 0.4) +
      stat_summary(alpha = 0.6,
                   fun   = mean,
                   geom  = "line") + # add a mean trendline per color


      # make sure the scale is the same on both axes, since both are in percent
      scale_y_continuous(breaks = seq(0, 1, 0.20), limits = c(0, 1)) +
      scale_x_continuous(breaks = seq(0, 100, 20)) +

      # add text
      ggtitle("Robustness of Method Predictions") +
      ylab("Overlap of Top Ranks [Jaccard Index]") +
      xlab(x_axis_name) +
      labs(subtitle = "Line / point scatter plot.",
           caption = plotting_caption,
           color = "Method") +

      # remove grey background so that points are bteer visible
      theme_bw() +

      # shift caption and legend position
      theme(plot.caption = element_text(hjust = 0),
            legend.position = "bottom",
            plot.margin = unit(c(0.5,0.5,1,0.5),"cm"))


    # return the plot
    return(plot_line)

  }

}


# overlap_box_Plot()
{

  #' Draw multiple boxplots plot from the collated top_ranks_overlap data
  #'
  #' @param top_ranks_overlap As a tibble, formatted like the collated_top_ranks
  #' overlap data.
  #'
  #' @param plotting_caption A caption for your plot, as a string. Potentially
  #' could be an auto_plot_Description() output.
  #'
  #' @return A ggplot of the top_ranks_overlap as a box plot, complete with a
  #' caption.

  overlap_box_Plot <- function(top_ranks_overlap,
                               plotting_caption,
                               x_axis_var,
                               x_axis_name) {

    plot_box <-
      ggplot(data = top_ranks_overlap, aes(x = .data[[x_axis_var]],
                                           y = Overlap,
                                           group = .data[[x_axis_var]], # for faceting
                                           color = Method)) +  # ...and colors

      # Add boxplots, no outliers because we add points in a second
      geom_boxplot(outlier.shape = NA) +
      geom_point(alpha = 0.4) +  # plot semi transparent plots over them so
      # you can see where the distirbution comes from


      # make sure the scale is the same on both axes, since both are in percent
      scale_y_continuous(breaks = seq(0, 1, 0.20), limits = c(0, 1)) +
      scale_x_continuous(breaks = seq(0, 100, 20)) +


      # add text
      ggtitle("Robustness of Method Predictions") +
      ylab("Overlap of Top Ranks [Jaccard Index]") +
      xlab(x_axis_name) +
      labs(subtitle = "Boxplot by Method.",
           caption = plotting_caption,
           color = "Method") +

      # remove gray background so the transparent points are more visible
      theme_bw() +

      # adjust legend and caption position
      theme(plot.caption = element_text(hjust = 0),
            legend.position = "bottom",
            plot.margin = unit(c(0.5,1,0.5,0.5),"cm")) +

      # facet wrap for each method, so they don't overlap and are easily
      # comparable between subplots
      facet_wrap(~Method, nrow = 3, ncol = 3, scales = "free")



    # return the plot we crafted
    return(plot_box)

  }


}



#------------------------------------------------------------------------------#
# 2. Define Functions for Console Outputs -------------------------------------

# print_Title()
{
  #' Prints a nice title to the console
  #'
  #' @param title A char of what the title should read
  #' @param super TRUE or FALSE. Should this title be extra large?
  #' @param space_after How many empty line should be after a header?
  #'
  #' @return Prints a title into the console.


  print_Title <- function(title, super = FALSE, space_after = 3) {
    divider <- str_glue(
      "|=======================================",
      "=======================================|"
    )

    if (super == FALSE) {
      cat(rep("\n", 3))
      print(str_glue(divider))
      cat(rep("\n", 1))
      print(str_glue("   ", title))
      cat(rep("\n", 1))
      print(str_glue(divider))
      cat(rep("\n", space_after))

    } else if (super == TRUE) {
      cat(rep("\n", 7))
      print(str_glue(divider))
      print(str_glue(divider))
      cat(rep("\n", 2))
      print(str_glue("   ", toupper(title)))
      cat(rep("\n", 2))
      print(str_glue(divider))
      print(str_glue(divider))
      cat(rep("\n", space_after))

    }



  } #end of function

}


# divert_Warnings()
{
  #' Replaces warnings from an expression with a logfile
  #'
  #' @description A function within which you can execute code. Warnings
  #' generated by the code are suppressed from the console and instead logged
  #' elsewhere.
  #'
  #' @param code Some expression to execute that produces warnings
  #' @param logFile The location you would like to log the warning at.
  #'
  #' @return A logfile with warnings instead of a console output with warnings.



  divert_Warnings <- function(code, logFile) {

    # Use a calling Handler to supress warning
    withCallingHandlers(code,
                        warning = function(w)
                        {
                          # cat warnings to a specified logFile
                          cat(as.character(Sys.time()), # time of warning
                              "\n",
                              as.character(w), # the warning as it would print
                              "\n",
                              file = logFile,
                              append = TRUE)

                          # don't print warning in console
                          invokeRestart("muffleWarning")
                        }
    )

  } # end of function

}

# 3. Reformat ----
#' Reformat plots to manuscript format
#' @param p plot to reformat
#'
#' @details plots come from the output of the robustness pipelines
#'
#' @returns a ggplot2 object
format_robustness_plot <- function(p, descript){
  p$data %<>%
    mutate(Method = ifelse(Method=="CytoTalk", "Crosstalk scores", Method))
  p +
    facet_grid(~Method, scales='free_x', space='free', switch="x") +
    theme_bw(base_size = 24) +
    theme(strip.text.x = element_text(angle = 90, face="bold", colour="white"),
          strip.background = element_rect(fill="darkgray"),
          legend.title = element_text(size = 28),
          legend.text = element_text(size = 25)
    ) +
    labs(colour=guide_legend(title="Method"),
         caption = NULL,
         title = NULL,
         subtitle = NULL) +
    ylab("Jaccard Index") +
    xlab(str_glue("{descript}")) +
    guides(shape = "none")
}



#' Function to plot TPR as Box plots
#'
#' @param res_to_tpr results converted to fpr/tpr (obtained via `convert_to_tpr`)
#' @param descirpt X-axis name
#'
#' @return a ggplot2 boxplot
plot_rates <- function(res_to_tpr, descript, rate='tpr'){
  if(rate=='tpr'){
    ylabel = "True Positive Rate"
  } else if(rate=='fpr'){
    ylabel = "False Positive Rate"
  }

  res_to_tpr %>%
    ggplot(aes(x = value,
               y = .data[[rate]],
               group = shuffle,
               colour = Method))  +
    geom_boxplot(outlier.shape = NA) +
    geom_point(alpha = 0.4) +  # plot semi transparent plots over them so
    # make sure the scale is the same on both axes, since both are in percent
    scale_y_continuous(breaks = seq(0, 1, 0.20), limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 100, 20)) +
    # add text
    ylab(ylabel) +
    labs(subtitle = "Boxplot by Method.",
         color = "Method") +
    facet_grid(~Method, scales='free_x', space='free', switch="x") +
    theme_bw(base_size = 24) +
    theme(strip.text.x = element_text(angle = 90, face="bold", colour="white"),
          strip.background = element_rect(fill="darkgray"),
          legend.title = element_text(size = 28),
          legend.text = element_text(size = 25)
    ) +
    labs(colour=guide_legend(title="Method"),
         caption = NULL,
         title = NULL,
         subtitle = NULL) +
    xlab(str_glue("{descript}")) +
    guides(shape = "none")

}

#' Function to calculate TPR from robustness results
#'
#' @param result_path path to the .rdata with the results from all runs.
#'
#' @return a tibble method, shuffle, and tpr
calculate_tpr <- function(result_path, analysis_focus){
  # load data
  load(result_path)

  # extract top ranks from all methods
  if(analysis_focus == "cluster"){
    top_ranks <- iterator_results$reshuffling_results$top_ranks
  } else if(analysis_focus == "resource"){
    top_ranks <- iterator_results$collated_robustness_results$top_ranks_OP
  }

  # calculate tpr in regards to 0 reshuffling
  res_to_tpr <- imap(top_ranks, function(method_res, method_name){

    # establish ground truth
    if(analysis_focus == "cluster"){
      gt_res <- top_ranks[[method_name]]$Reshuffle_0$Seed_1
    } else if(analysis_focus == "resource"){
      gt_res <- top_ranks[[method_name]]$OmniPath_0_Seed_1
    }

    method_res %>%
      enframe(name='shuffle', value = 'results') %>%
      {`if`(analysis_focus == 'cluster', unnest(., results), .)} %>%
      mutate(jacc = map_dbl(results, function(res){
        rank_overlap(comparison_ranks = res, main_ranks = gt_res)
      })
      ) %>%
      mutate(tpr = map_dbl(results, function(res){
        tpr_overlap(comparison_ranks = res, main_ranks = gt_res)
      })
      ) %>%
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

  return(res_to_tpr)
}

