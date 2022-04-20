# Install renv
install.packages("renv")

# Use renv to restore environment
renv::restore()

# Install the appropriate liana version
remotes::install_github("saezlab/liana", ref="fdb2018")

# Install LIANA++ version following this tutorial:
# https://saezlab.github.io/liana/articles/liana_devel.html

# Obtain OmniPath /w mouse symbols and save
source("src/eval_utils.R")
murine_resource <- liana::select_resource("OmniPath")[[1]] %>%
    convert_to_murine()
saveRDS(murine_resource, "data/input/murine_omnipath.RDS")

