# Bash Scripts Used to Generate LIANA results
This directory contains the bash and slurm job scripts used to orchestrate the submissions of
tasks which require long run times (e.g. the comparison of all resource-method combinations)
or such which require a large amount of RAM.

1) Use `main_bash.sh` after all Seurat objects were prepared to launch all LIANA runs.

2) Use `comp_summarise.sh` to launch all comparison summary scripts

3) Use `extract_evals.sh` to extract and aggregate OmniPath results for each method for all BRCA subtypes

Note: generic- large- and huge_bash are identical scripts with different RAM and wallclock requirements.
