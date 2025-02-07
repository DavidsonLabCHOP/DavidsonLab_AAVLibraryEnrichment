#!/bin/bash

####################################
## This wrapper is designed to be run from a folder that contains both
##    your data files and the AAV_Barcode_Parsing_Tool.py script.

########################################
## Necessary User-Supplied Input
## 
## You must supply a text file containing the names of the experimental
##    data files (*.fastq) that you want to analyze. This sample wrapper script
##    is written assuming this text file is in the same folder as this script.
## Here we've called our text file FASTQ_Filenames.txt 
## If you save your list of experimental data files to a text file with a different
##    name, change the script below.

mapfile -t key_Array < FASTQ_Filenames.txt

for i in "${key_Array[@]}"
do
  SampleName=$(echo $i | cut -d'.' -f1)
  echo "Sample Name is: $SampleName"

  echo "beginning counts generation workflow"
  python3 AAV_Barcode_Parsing_Tool.py -D "." -F "$i" -A true > "$SampleName-counts.txt"
  
done
