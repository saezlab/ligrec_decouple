# Systematic Comparison of Cell-Cell Communication Tools and Resources

Repository used to reproduce the results from [Comparison of Resources and Methods to infer Cell-Cell Communication from Single-cell RNA Data](https://www.biorxiv.org/content/10.1101/2021.05.21.445160v1)


## LIANA Analysis Content

### I) Descriptive Resource Analysis
The code to reproduce the descriptive analysis of resources can be found at:
[analysis/resource_analysis](https://github.com/saezlab/ligrec_decouple/tree/main/analysis/resource_analysis)

### II) Spatial Co-localization
The code to reproduce the co-localization analysis can be found at:
[analysis/spatial](https://github.com/saezlab/ligrec_decouple/tree/main/analysis/spatial)

### III) Cytokine Signalling Agreement
The code to reproduce the cytokine activity (/w [CytoSig](https://www.nature.com/articles/s41592-021-01274-5)) agreement analysis can be found at:
[analysis/cytosig](https://github.com/saezlab/ligrec_decouple/tree/main/analysis/cytosig)

### IV) CITE-Seq Correlation/Specificity
The code to reproduce the Correlation/Specificity analysis of methods with CITE-Seq can be found at:
[analysis/citeseq](https://github.com/saezlab/ligrec_decouple/tree/main/analysis/citeseq)

### V) Comparison of Methods and Resources
The code to reproduce the comparison between method-resource combinations can be found at:
[analysis/comparison](https://github.com/saezlab/ligrec_decouple/tree/main/analysis/comparison)

### VI) Robustness
The code to reproduce the robustness analyses can be found at:
[analysis/robustness](https://github.com/saezlab/ligrec_decouple/tree/main/analysis/robustness)


## Environment set-up
# Clone repo
```{bash}
git clone https://github.com/saezlab/ligrec_decouple
```

Use `initiate_environment.R` to restore the R environment used to produce the results in the publication.

Finally, refer to [LIANA++](https://saezlab.github.io/liana/articles/liana_devel.html) for instructions how to set-up the LIANA environment used to call the different external tools.

