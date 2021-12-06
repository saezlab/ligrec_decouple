## Submit Robustness Runs

### Cluster Reshuffling
bash bash_scripts/large_bash.sh "analysis/robustness/Code/Analysis_Scripts/Clusters_Run_Iterator.R" "clustre" "reshuffle" $SLURM_JOBID

### Cluster Subsetting
bash bash_scripts/large_bash.sh "analysis/robustness/Code/Analysis_Scripts/Clusters_Run_Iterator.R" "clussub" "subset" $SLURM_JOBID

### Indiscriminant Resource dilution (Modify baseline = TRUE)
bash bash_scripts/large_bash.sh "analysis/robustness/Code/Analysis_Scripts/Resources_Run_Iterator.R" "indscres" TRUE $SLURM_JOBID

### Discriminant Resource dilution (Modify baseline = FALSE)
bash bash_scripts/large_bash.sh "analysis/robustness/Code/Analysis_Scripts/Resources_Run_Iterator.R" "discres" FALSE $SLURM_JOBID
