#!/bin/bash
## declare an array variable
declare -a subtypes=("ER2" "TNBC" "HER2")

# get length of an array
arraylength=${#subtypes[@]}

# use for loop to read all values and indexes
for (( i=0; i<${arraylength}; i++ ));
do
 Rscript ./brca_prep.R ${subtypes[$i]}
done
