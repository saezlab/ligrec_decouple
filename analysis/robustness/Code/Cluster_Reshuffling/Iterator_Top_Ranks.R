#------------------------------------------------------------------------------#
# 0. Introduction and Goals ----------------------------------------------------
{
  # This is the Iterator_Top_Ranks.R script. It defines functions the Iterator
  # uses to manipulate top-ranked CCIs.
}



#------------------------------------------------------------------------------#
# 1. Define Functions ----------------------------------------------------------

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



