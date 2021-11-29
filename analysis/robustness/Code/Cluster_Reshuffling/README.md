---
editor_options: 
  markdown: 
    wrap: 72
---

# Cluster_Reshuffling:

Welcome! This is the code folder for the cluster reshuffling process.
The goal of this script is to measure the robustness of CCI (Cell-Cell
Interaction) inference methods when the cluster annotations of your
data become more and more inaccurate (or sparse).

The most important script is the **Run_Iterator.R** script. It executes
a large wrapper that performs the entire analysis. All the other code
simply defines the wrapper or subfunctions used in the wrapper. So, if
you just want to run an analysis, go to **Run_Iterator.R**. It also has
a more detailed overview of the theoretical approach used in the
analysis.

If you want more insight, feel free to check the other scripts too.

# Script TOC:

1.  **Run_Iterator.R**: This script that executes the wrapper for the
    analysis.

2.  **CR_Iterator.R**: This script defines the wrapper for the analysis.

3.  **Iterator_Top_Ranks.R**: This script defines helpers for the
    iterator.

4.  **Iterator_Meta_and_Saves.R**: This script defines helpers for the
    iterator.

5.  **CR_Shuffler_Functions.R**: These functions will reshuffle cluster
    annotations.

6.  **CR_Subsetter_Functions.R**: These functions will subset cluster
    annotations.

7.  **CR_LIANA_Functions.R**: These functions run LIANA over reshuffled
    or subset cluster annotations.
