#!/bin/bash
## The one bash script to control them all

## Submit all BRCA scripts
bash comp_brca_prep ER
bash comp_brca_prep TNBC
bash comp_brca_prep HER2

## Submit Comparisons

## Submit CITE-Seq pipe
