# DavidsonLab_AAVLibraryEnrichment
Tools for processing AAV screening data from Illumina amplicon sequencing

This repository contains a series of bash/python/R tools for analyzing Illumina next-generation amplicon sequencing results from peptide-modified AAV screening experiments.
The code is designed to analyze data from amplicons that contain heptamer (21 base pair) variable regions. Properly formed amplicons in the Illumina FASTQ results will identified using perfect string matches. This requires that the variable region, two flanking sequences, and a Unique Molecular Identifier sequence (UMI) be at fixed locations in the amplicon. All of these sequences must be located on the same Read in the Illumina data.

---------------------
AMPLICON FASTQ PARSER
---------------------
This folder contains the python script that will parse the Illumina FASTQs: [AAV_Barcode_Parsing_Tool.py]. And an sample bash wrapper script for calling the python tool: [AAV_Parsing_Tool_Wrapper.sh].

Annotations, including important instructions for running these scripts, can be found in the scripts themselves. Below are the minimum requirements necessary before running these scripts.

AAV_Parsing_Tool_Wrapper.sh
---------------------------
Before running this wrapper, users will need to create a text file containing a list of the FASTQs they want to analyze. Users will then need to modify the AAV_Parsing_Tool_Wrapper.sh script to replace the <FASTQ_Filenames.txt> placeholder with the name of their own text file.
If no further modifications are made to the AAV_Parsing_Tool_Wrapper.sh script, then all of the following must be located in the same folder prior to running the wrapper: 
1) AAV_Parsing_Tool_Wrapper.sh
2) The AAV_Barcode_Parsing_Tool.py script
3) The FASTQs to be analyzed (FASTQs must be unzipped)
4) The text file containing a list of filenames for the FASTQs to be analyzed

AAV_Barcode_Parsing_Tool.py
---------------------------
This script is designed to be run from a Linux/UNIX Shell script wrapper. Before running this script, there are 9 variables that require user customization. All of these variables can be found at the top of the script:
1) flank_1_start_index - The position of the first nucleic acid in the Flank1 sequence
2) flank_1_end_index - The position of the last nucleic acid in the Flank1 sequence
3) flank_1_sequence - A string defining the expected nucleic acid sequence in Flank1
4) flank_2_start_index - The position of the first nucleic acid in the Flank2 sequence
5) flank_2_end_index - The position of the last nucleic acid in the Flank2 sequence
6) flank_2_sequence - A string defining the expected nucleic acid sequence in Flank2
7) peptide_modification_start_index - The position of the first nucleic acid in the heptamer peptide modification
8) UMI_start_index - The position of the first nucleic acid in the UMI
9) UMI_end_index - The position of the last nucleic acid in the UMI

This script will output a <original_filename>-counts.txt file for each Illumina FASTQ analyzed. These counts files will contain raw read counts and UMI collapsed read counts for each unique variable region sequence detected in the Illumina FASTQ

---------------------------------------
VARIANT_PERFORMANCE_ANALYSIS
---------------------------------------
This folder contains an R script that analyzes the *counts.txt files output from the Amplicon_FASTQ_Parser tools to evaluate capsid variant performance using two metrics: UMI counts and Enrichment.

AAV_Screening_Data_Analysis_Tool.R
----------------------------------
Complete instructions for running this tool along with explanations of all available options and outputs can be found in the annotations contained within the script.

The only input required to run this analysis tool is a text file (UTF-8 encoded) containing the list of *-counts.txt files to be analyzed. This text file must be located in the same folder as the *-counts.txt files. Also note, one of the *-counts.txt files must contain data for the input virus. Input virus data is required to calculate the Enrichment metric.







