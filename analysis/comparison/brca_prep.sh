#!/bin/bash
## declare an array variable
declare -a subtypes=("ER" "TNBC" "HER2")

# get length of an array
arraylength=${#subtypes[@]}

# Loop over brca_types and first comm arg is path_to_project
for (( i=0; i<${arraylength}; i++ ));
do
 echo "Submitted: Rscript ./comp_brca_prep.R $1 ${subtypes[$i]}"
 Rscript ./comp_brca_prep.R $1 ${subtypes[$i]}
done
