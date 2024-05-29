#!/bin/bash

# first accept a path to the json parameter files of interest
read -p 'Enter path to directory with json param files: ' param_dir
for filename in "$param_dir"/*; do
  echo "Now simulating with parameter file: $filename"
  Rscript sim5_code.R -P "$filename"

done
