# Install renv
install.packages("renv")

# Use renv to restore environment
renv::restore()

# Obtain OmniPath /w mouse symbols and save
source("src/eval_utils.R")
murine_resource <- liana::select_resource("OmniPath")[[1]] %>%
    convert_to_murine()
saveRDS(murine_resource, "data/input/murine_omnipath.RDS")
