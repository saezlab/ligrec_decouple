## LIANA Analysis Content

Use `initiate_environment.R` to restore the R environment used to produce the results in the publication.
Refer to [LIANA++](https://saezlab.github.io/liana/articles/liana_devel.html) for instructions how to set-up the LIANA environment used to call the different external tools.

### I) Descriptive Resource Analysis
The code to reproduce the descriptive analysis of resources can be found at:
[analysis/resource_analysis](https://github.com/saezlab/ligrec_decouple/tree/main/analysis/resource_analysis)

### II) Comparison of Methods and Resources
The code to reproduce the comparison between method-resource combinations can be found at:
[analysis/comparison](https://github.com/saezlab/ligrec_decouple/tree/main/analysis/comparison)

Note that some of the other types of analyses make use of the output of the comparisons.
Thus, it should be ran before the others are.

### III) Spatial Co-localization
The code to reproduce the co-localization analysis can be found at:
[analysis/spatial](https://github.com/saezlab/ligrec_decouple/tree/main/analysis/spatial)

### IV) Cytokine Signalling Agreement
The code to reproduce the cytokine activity (/w [CytoSig]()) agreement analysis can be found at:
[analysis/cytosig](https://github.com/saezlab/ligrec_decouple/tree/main/analysis/cytosig)

### V) CITE-Seq Correlation/Specificity
The code to reproduce the Correlation/Specificity analysis of methods with CITE-Seq can be found at:
[analysis/citeseq](https://github.com/saezlab/ligrec_decouple/tree/main/analysis/citeseq)

### VI) Robustness
The code to reproduce the robustness analyses can be found at:
[analysis/robustness](https://github.com/saezlab/ligrec_decouple/tree/main/analysis/robustness)


Note that `plot_evals_figure6.R` is the script ran on both CytoSig and Spatial output to 
produce the figure presented in the manuscript.
