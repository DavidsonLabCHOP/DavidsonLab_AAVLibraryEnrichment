# This script is set up to process AAV Evolution Data - Round 2 Data for this build
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(scales)
library(officer)
library(rvg)
library(purrr)
library(sys)
library(reshape2)

##############################################
##  Important code settings and variables   ##
##############################################

###############################################
###  Data Processing and Normalization Options  ###
###############################################

## This toggle will cause the script to delete all the variables it creates once it has finished
#     running and free up any unused memory. 
Clear_Environment_Upon_Completion <- TRUE

## This variable defines the UMI count threshold used to filter out weakly detected Barcodes.
#     This value indicates the minimum number of UMIs required for a given Barcode to be included 
#     in enrichment analyses.
UMI_Per_Barcode_Threshold <- 10





#######################################
##
##  Sequencing Error Detection and Removal Options

## This toggle determines whether the script will perform automated sequence error barcode detection
##    and removal. Note: This is a memory intensive process. It is run individually on each dataset. 
##    If a given dataset includes too many barcodes, it is possible that R could run out of available memory. 
##    The UMI_Per_Barcode_Threshold is applied before the sequence error detection/removal is run. So, you can
##    try increasing that threshold if you are having memory issues.
Collapse_Barcodes <- TRUE

## This variable determines the edit distance between barcodes used to determine which barcodes
##    are compared as potential sequencing errors.
edit_dist_setting <- 1

## Putative sequencing error barcodes need to have less than this Percent_of_Parent threshold 
#     in order to be considered an actual sequencing error barcode. Percent of Parent is calculated
#     using the following equation for each putative sequencing error barcode
#     (Err_Barcode_UMIs / Parent_Barcode_UMIs) * 100
#     Note: this value is a % (i.e., it ranges from 0 to 100, not 0 to 1).
Percent_of_Parent <- 1





## The string saved to the Save_File_Tag variable will be appended to the end of the standard names for each of the outputs
#     saved by this script. You'll want to make this unique to each data set so you don't accidentally over-write
#     the output of a previous run.

Save_File_Tag = "Your_Experiment_ID"







#############
############################
##################################################3
#############################
#############

Data_File <- "Full_Path_to_Text_File_Containing_List_of_Counts_Files_to_be_Analyzed"

  
Data_File_List <- basename(Data_File)
Experimental_Data_Dir <-  dirname(Data_File)
setwd(Experimental_Data_Dir)
  
  







#######################################
###  Identification of Input Vectors  ###
#######################################

### Normalization and data analysis in this script requires that the input vector files be analyzed separately 
#     from the experimental results files. You can enter the Input Vector file IDs manually or try an automated
#     identification of input files by setting AutoFind_Input_Files = TRUE.
#     The auto-find feature relies on the Input Vector files having "Input" in their file name.

AutoFind_Input_Files <- TRUE

# You only need to define the Input_Vector_IDs below if the AutoFind_Input_Files toggle is set to FALSE
#AAV_Input_Vector_ID <- 'Carrells-2-NA-Vector-R2-Input-AAV1-NA-Universal-Merged-v1'




#############################
###  Script Output Options  ###
#############################

## This option will save a list structure of results tables for each dataset. These tables includes the 
##   barcode, peptide insert, normalized experimental data, performance metric and normalized input vector data.
Save_Results_Tables = TRUE

## Change this to choose where the results are saved. Use getwd() to save in same location as the experimental files.
##    Currently, the script creates a folder called Results in the folder where the counts files are located to save the output files
Save_Dir = paste0(getwd(),"/Results","-",Save_File_Tag)
dir.create(Save_Dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



###############################
###   STEP1: Import Data    ###
###############################

# These first five lines of code set up a couple of variables that will track certain activity in the
#   script (e.g., removal of problematic datasets). The information they store is output as part of the
#   Settings_and_Info text file that is saved at the end of the script.
QC_Report <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(QC_Report) <- c("Data_File","Reason_For_Exclusion")

Output_Notes <<- array()
Output_Notes[1] <- "Additional Information:"
Output_Notes[2] <- "-----------------------"


# This next set of code is going to import the data files created by the Python scripts run on the amplicon FASTQ files

# These datasets contain four columns of data:
# Column 1 is the Barcode sequence in DNA base pairs
# Column 2 is the Barcode sequence translated to amino acid letter code
# Column 3 is the total counts
# Column 4 is the UMI normalized counts


### Create a list of sample IDs
SampleIDs <- data.frame(read_csv(Data_File_List, col_names=FALSE))
names(SampleIDs) <- c("Names") 


##################################################################
#### Create a function to import the dataset of interest        ##
##################################################################

importdatasets <- function(sample) {
  outputDF <- data.frame(read_table(sample, col_names = FALSE ))
  return(outputDF)
}

### Apply the importdatasets function to read in the datasets defined by the ID_list.txt file
datasets=list()
datasets <- lapply(SampleIDs$Names, FUN = importdatasets)
names(datasets) <- SampleIDs$Names  

for (i in c(1:length(datasets))) {
  colnames(datasets[[i]]) <- c("Barcodes","Peptide","RawCounts","UMICounts")
}


# Calculate the percent of total reads in a dataset accounted for by each barcode
alldatasets <- lapply(datasets, function(x) {x["pctTotalReads"] <- 100 * (as.numeric(x$RawCounts) / sum(as.numeric(x$RawCounts))); return(x)})






