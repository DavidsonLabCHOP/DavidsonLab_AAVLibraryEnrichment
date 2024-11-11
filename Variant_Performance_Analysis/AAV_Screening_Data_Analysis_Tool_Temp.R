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
UMI_Per_Barcode_Threshold <- 100





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
##################################################
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
#AAV_Input_Vector_ID <- 'Name of your input vector counts.txt file'




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
alldatasets=list()
alldatasets <- lapply(SampleIDs$Names, FUN = importdatasets)
names(alldatasets) <- SampleIDs$Names  

for (i in c(1:length(alldatasets))) {
  colnames(alldatasets[[i]]) <- c("Barcodes","Peptide","RawCounts","UMICounts")
}


# Calculate the percent of total reads in a dataset accounted for by each barcode
alldatasets <- lapply(alldatasets, function(x) {x["pctTotalReads"] <- 100 * (as.numeric(x$RawCounts) / sum(as.numeric(x$RawCounts))); return(x)})











#######################################################
###   AAV Dataset Processing                        ###
#######################################################



###############
#############################
#  Identify and Separate Input Vector Data from List of Experimental Data Files
#
#  To improve signal-to-noise in capsid performance evaluation, and as a first step in removal of potential sequencing errors,
#     datasets are filtered to remove barcodes with low detection levels. 
#  Unlike the experimental datasets, we want to preserve all information from the input virus. Its data will be used to
#     calculate enrichment metrics for barcodes present in experimental datasets.

if (AutoFind_Input_Files) {
  AAV_Input_Vector_ID <- names(alldatasets[grep("Input", names(alldatasets))])
}
input_vector_dataset <- alldatasets[[AAV_Input_Vector_ID]]
alldatasets[[AAV_Input_Vector_ID]] <- NULL

# Adding input vector file ID info to the run-info variable
Output_Notes[length(Output_Notes)+1] <- paste("AAV Input Vector File ID:",as.character(AAV_Input_Vector_ID)) 



###############
#############################
#  Apply UMI_Per_Barcode Threshold filter

## Filter datasets to remove barcode entries with UMI counts < UMI_Per_Barcode_Threshold
alldatasets <- lapply(alldatasets, function(x) { subset(x, x$UMICounts >= UMI_Per_Barcode_Threshold) })

## Identify, track and remove datasets with no remaining barcodes following UMI thresholding
for (i in length(alldatasets):1) {
  if (nrow(alldatasets[[i]])==0) {
    ## The QC_Report data.frame is used to track datasets that are excluded and save that information to the 
    #       Settings_And_Info.txt file output at the end of the script
    QC_Report[nrow(QC_Report)+1,] <- c(names(alldatasets)[i],"Dataset Empty after UMI Filter Applied")
    alldatasets <- alldatasets[-i]
  }
}



###############
#######################################
#  Collapse Barcodes
#
#  This function collapses barcodes based on amino acid sequence distance. The toggle that controls
#    whether to collapse barcodes and the variable that controls the edit distance threshold are
#    located at the top of the script.

if (Collapse_Barcodes) {
  
  ## Note: We don't want to perform sequence error removal on the input vector, so we aren't adding it back to the 
  #     dataset list yet. 
  #  Only barcodes present in experimental data files will be analyzed; so, it is a waste of processing time to  
  #     perform sequencing error removal on the input vector.
  #  In addition, it is likely the input vector contains enough barcodes that attempting to process it could
  #     use all available memory and crash R.
  
  alldatasets <- lapply(alldatasets, collapse_BCs)
}


## Adding the input virus data back to the list of experimental datasets
alldatasets[[AAV_Input_Vector_ID]] <- input_vector_dataset






########################################################################################
##   Create a master key for all barcode sequences detected across all datasets     ####
########################################################################################

# This code is going to create one monster data frame that combines all the individual data
#     frames, stacking them one on top of the other. The goal is to create a reference for
#     every barcode and peptide sequence observed across all the datasets for a given serotype.
#     We can then use this "Barcode/Peptide" key to retrieve this information later on
#     after we have processed the data in ways that could filter out some of this info.

AAV_Final_Key <- do.call("rbind", alldatasets)
AAV_Final_Key <- AAV_Final_Key[,c('Barcodes','Peptide')]

#   Now we're going to remove duplicate barcodes and make peptide sequences unique
AAV_Final_Key <- AAV_Final_Key[!duplicated(AAV_Final_Key$Barcodes),]
AAV_Final_Key["Peptide"] <- make.unique(AAV_Final_Key$Peptide)









######################################################################################
######################################################################################
# 
#   IMPLEMENT NORMALIZATION METHOD
#   
#   The code in this section normalizes the data. There are multiple options
#   for how the data can be normalized. In its current incarnation, this code
#   runs one normalization function and also passes through unnormalized data for 
#   some control calculations. Normalization is implemented by calling one of the 
#   NormFun functions. This call is made from within the NormalizationAndFirstPassRankings function
#
#   The NormalizationAndFirstPassRankings function is called by the code found in the section after
#   this one.
#
#######################################################################################
#######################################################################################



Create_Metric_Data_Frames <- function(input_vector, experimental_data) {
  
  ## Merging the experimental and input virus datasets
  #    Note - This merging function removes barcodes from the experimental dataset that are not detected in the
  #      input vector. You will want to run QC checks to ensure that this process is not removing any barcodes of
  #      potential interest. If it is, you may need to improve the sequencing depth of your input vector library.
  experimental_data <- as.data.frame(merge(input_vector[,c('Barcodes','Peptide','UMICounts')], experimental_data[,c('Barcodes','UMICounts')], by = "Barcodes"))
  colnames(experimental_data) <- c("Barcodes", "Peptide", "InputVector", "ExpData")
  
  experimental_data[is.na(experimental_data)] <- 0
  
  
  experimental_data_enrichment <- experimental_data[,c('Barcodes','Peptide')]
  experimental_data_enrichment$PctUMI_Input <- experimental_data$InputVector / sum(experimental_data$InputVector)
  experimental_data_enrichment$PctUMI_Exp <- experimental_data$ExpData / sum(experimental_data$ExpData)
  experimental_data_enrichment$Enrichment <- (experimental_data_enrichment$PctUMI_Exp - experimental_data_enrichment$PctUMI_Input)/experimental_data_enrichment$PctUMI_Input
  experimental_data_enrichment$Enrichment <- experimental_data_enrichment$Enrichment + 1
  
  experimental_data_enrichment <- dplyr::arrange(experimental_data_enrichment, desc(Enrichment))
  experimental_data_enrichment$Enrichment_Rank <- c(1:nrow(experimental_data_enrichment))
  
  
  colnames(experimental_data) <- c("Barcodes","Peptide","UMI_Input","UMI_Exp")
  experimental_data$UMI_Input <- as.numeric(experimental_data$UMI_Input)
  experimental_data$UMI_Exp <- as.numeric(experimental_data$UMI_Exp)
  
  experimental_data <- dplyr::arrange(experimental_data, desc(UMI_Exp))
  experimental_data$UMIRank_Exp <- c(1:nrow(experimental_data))
  
  return(list(experimental_data_enrichment, experimental_data))
}





#########################################################################################
#########################################################################################
##
##    Calculate Enrichment Metric Values - Save Basic Metric Data Frames
##
#########################################################################



###################################
# Run Datasets                    #
###################################

setwd(Save_Dir)

# This variable is used in the Settings_And_Info output file
saved_files <- array()
saved_files[1] <- "Names of Saved Datafiles:"
saved_files[2] <- "----------------------"

##  This counter is used to provide some minor progress update functionality
counter <- 1

## We've already removed datasets with no data in them. However, the Create_Metric_Data_Frames function will 
#     crash if there are datasets with only a single row of data. This code identifies and removes datasets with
#     only a single barcode. Incidentally, a dataset having only 1 barcode is almost certianly an indication
#     of very low quality data.

for (i in length(alldatasets):1) {
  if (nrow(alldatasets[[i]])<=1) {
    QC_Report[nrow(QC_Report)+1,] <- c(names(alldatasets[i]),"Only 1 Barcode in Dataset")
    alldatasets <- alldatasets[-i]
  }
}



## The code below creates a simple function whose primary purpose is to call the NormalizationAndFirstPassRankings function. 
#     We then call this function using lapply to run all the datasets through the normalization and analysis.
#     There is a print command linked to a counter that will display the file being processed. This could be commented out
#     to speed things up a bit.

RunAAV_Samples <- function(exp_data) {
  print(paste0("Processing ",as.character(names(alldatasets)[counter])))
  processed_data <- Create_Metric_Data_Frames(input_vector_dataset, 
                                              exp_data)
  
  counter <<- counter + 1
  return(processed_data)
}


Metric_Tables <- lapply(alldatasets, RunAAV_Samples)



## Save the enrichment and UMI counts data frames created by RunAAV_Samples
if (Save_Results_Tables) {
  
  Extract_Enrichment_Tables <- function(EnrichmentData) {
    EnrichmentData <- data.frame(EnrichmentData[[1]])
    return(EnrichmentData)
  }
  Extract_UMIcount_Tables <- function(UMIData) {
    UMIData <- data.frame(UMIData[[2]])
    return(UMIData)
  }
  
  Enrichment_Tables <- lapply(Metric_Tables, Extract_Enrichment_Tables)
  UMIcount_Tables <- lapply(Metric_Tables, Extract_UMIcount_Tables)
  
  saveRDS(Enrichment_Tables, file = paste0(as.character(Save_Dir),"/Enrichment_Tables-",Save_File_Tag,".rds"))
  saveRDS(UMIcount_Tables, file = paste0(as.character(Save_Dir),"/UMIcount_Tables-",Save_File_Tag,".rds"))
  saveRDS(AAV_Final_Key, file = paste0(as.character(Save_Dir),"/AAV_Final_Key-",Save_File_Tag,".rds"))
  
  saved_files[length(saved_files)+1] <- paste0("Enrichment_Tables-",Save_File_Tag,".rds")
  saved_files[length(saved_files)+1] <- paste0("UMIcount_Tables-",Save_File_Tag,".rds")
  saved_files[length(saved_files)+1] <- paste0("AAV_Final_Key-",Save_File_Tag,".rds")
  
}





General_Output_Filename <- paste0("Settings_And_Info-",Save_File_Tag,".txt")

sink(file = General_Output_Filename)

print("AAV_Screening_Data_Analysis_Tool.R - Setting and Information")
print("Run on:")
Sys.time()
print("---------------------------------------------------------------------------------------------")
print("---------------------------------------------------------------------------------------------")
print("Files and Directories")
print("---------------------")
print(paste("Data File List:",as.character(Data_File_List)))
print(paste("Experimental Data Directory:",as.character(Experimental_Data_Dir)))
print(paste("Save Directory:", as.character(Save_Dir)))
print(paste("Save File Identifier:",Save_File_Tag))
print("---------------------------------------------------------------------------------------------")
print("---------------------------------------------------------------------------------------------")
print("Normalization and Filtering Parameters")
print("--------------------------------------")
print(paste("UMI_Per_Barcode_Threshold =",as.character(UMI_Per_Barcode_Threshold)))
print(paste("Collapse_Barcodes =",Collapse_Barcodes))
print(paste("edit_dist_setting =",as.character(edit_dist_setting)))
print(paste("consolidation_percentage =",as.character(consolidation_percentage)))
print("---------------------------------------------------------------------------------------------")
print("---------------------------------------------------------------------------------------------")
print("Option Toggles")
print("--------------")
print(paste("AutoFind_Input_Files =",AutoFind_Input_Files))
if (!AutoFind_Input_Files) {
  print(paste("AAV_Input_Vector_ID:",AAV_Input_Vector_ID))
}
print(paste("Save_Results_Tables =", Save_Results_Tables))
print("---------------------------------------------------------------------------------------------")
print("---------------------------------------------------------------------------------------------")
print("Files Analyzed:")
print("-------------")
for (i in 1:nrow(SampleIDs)) {
  print(SampleIDs$Names[i])
}
print("---------------------------------------------------------------------------------------------")
print("---------------------------------------------------------------------------------------------")
for (i in 1:length(saved_files)) {
  print(saved_files[i])  
}
print("---------------------------------------------------------------------------------------------")
print("---------------------------------------------------------------------------------------------")
print("QC Report File")
print("--------------")
QC_Report
print("---------------------------------------------------------------------------------------------")
print("---------------------------------------------------------------------------------------------")
print("Output Notes:")
print("-------------")
for (i in 1:length(Output_Notes)) {
  print(Output_Notes[i])  
}
sink(file = NULL)


General_Output_Filename <- paste0(Save_File_Tag,"-AllFiles.txt")

sink(file = General_Output_Filename)
for (i in 3:length(saved_files)) {
  print(saved_files[i])  
}
sink(file = NULL)




























#########################################################################################
#########################################################################################
##
##    COLLAPSE BARCODE FUNCTION
##
##    The toggle to control whether this function is run can be found at the beginning
##        of the script. In addition, the variable that determine the error rate (edit_dist_setting)
##        can also be found at the top of the script.
##
##    This function will find barcodes within edit_dist_setting of each other and combine
##        then into the a single entry based on which one of them has the highest UMI counts.
##
#########################################################################


collapse_BCs <- function(data) {
  
  ## I'm renaming the column rows here because I used slightly different names when
  #     initially writing this function.
  colnames(data) <- c('Barcodes','Peptide','RawCounts','UMIs','pctTotalReads')
  
  ###############################################
  ##  The first steps are the same as before. Calculate edit_distances and then identify
  ##      the oPool barcode pairs within 'edit_dist_setting' of each other.
  
  # Calculate the pairwise edit distance between strings
  # If you are looking for the step that takes the longest to run in this code, this is it. The adist()
  #   function creates a pair-wise matrix which grows exponentially with the number of comparisons (Barcodes).
  # This adist() function takes a list variable to determine the costs of various string manipulations.
  #   We are only interested in substitutions. To prevent insertion or deletion based results registering as
  #   valid barcode comparisons, the cost of insertions and deletions is set to 10.
  #       In contrast, the cost of substitutions is set to 1, so the value of the edit_dist_setting variable  
  #       will be directly correlated with the maximum allowable number of edits.
  edit_costs <- list()
  edit_costs[['insertions']] <- 10
  edit_costs[['deletions']] <- 10
  edit_costs[['substitutions']] <- 1
  
  edit_distances <- adist(data$Barcodes,costs=edit_costs,counts=TRUE)
  
  # Filter the pairs with an edit distance of <edit_dist_setting> or less (excluding the diagonal elements)
  pairs <- which(edit_distances <= edit_dist_setting & edit_distances > 0, arr.ind = TRUE)
  
  
  #################################################
  ##  Creating a new data frame, barcode_similarities, to store the count of other
  #     barcodes within <edit_dist_setting> of each individual barcode. This data frame is just a list
  #     of the barcodes in 'data' and a column called 'BCs_within_edit_dist' that contains that info.
  
  ## Create a new empty data frame and populate it with the barcode info from 'data'
  barcode_similarities <- data.frame(matrix(nrow=nrow(data),ncol=2))
  colnames(barcode_similarities) <- c('Barcodes','BCs_within_edit_dist')
  barcode_similarities$Barcodes <- data$Barcodes
  
  for (i in 1:nrow(data)) {
    ## Pull information on other barcodes within <edit_distance_setting> of the current barcode from the
    #     pairs matrix.
    barcode_similarities$BCs_within_edit_dist[i] <- length(which(pairs[,2] == i))
  }
  
  ## I wanted to be able to see the BCs_within_edit_dist data as part of the output of the
  #     following functions, so I added that info to the 'data' variable. That way I can
  #     use this new data_with_EdDistCnt data frame with rbind in the functions below and get
  #     outputs that have all the 'data' info plus the BCs_within_edit_dist.
  data_with_EdDistCnt <- cbind(data,barcode_similarities$BCs_within_edit_dist)
  colnames(data_with_EdDistCnt) <- c('Barcodes','Peptide','RawCounts','UMIs','pctTotalReads','BCs_within_edit_dist')
  #data_with_EdDistCnt <- data_with_EdDistCnt[order(data_with_EdDistCnt$BCs_within_edit_dist,decreasing = TRUE),]
  
  #####################################################################################
  ##
  ##  Create a data frame containing the "parent" barcodes that should be used as the 
  ##      primary input that the other barcodes will be collapsed into
  ##
  ##  This function identifies "parent" barcodes based on which barcode has the
  ##      highest UMI counts in pairwise comparisons.
  ##
  ######################################################################################
  
  #Test_BCcons_UMI_Dist <- data.frame(matrix(nrow=0,ncol=2))
  #colnames(Test_BCcons_UMI_Dist) <- c('Barcodes','Avg_UMI_Diff')
  
  BCcons_Indiv_Comps <- list()
  cross_comparisons <- list()
  
  BCs_to_consolidate_UMI <- data.frame(matrix(nrow=0,ncol=6))
  colnames(BCs_to_consolidate_UMI) <- c('Barcodes','Peptide','RawCounts','UMIs','pctTotalReads','BCs_within_edit_dist')
  
  for (i in 1:nrow(data)) {
    
    Test_BC_UMIs <- data_with_EdDistCnt$UMIs[i]
    
    ## We still need to grab the BCs_within_edit_dist value for our test barcode because
    #     we don't want to waste time checking data for unique barcode entries.
    Test_BC_EdDistCount <- data_with_EdDistCnt$BCs_within_edit_dist[i]
    
    if (Test_BC_EdDistCount > 0) {
      
      bcs_within_edit <- which(pairs[,2] == i)
      bcs_within_edit <- pairs[bcs_within_edit,1]
      
      ## This is a little different than the previous function because we need to grab the UMI
      #     data instead of the BCs_within_edit_dist data.
      cross_comparison <- data.frame(Barcodes = matrix(nrow=length(bcs_within_edit),ncol=1))
      cross_comparison$Barcodes <- barcode_similarities$Barcodes[bcs_within_edit]
      cross_comparison <- left_join(cross_comparison,data_with_EdDistCnt[,c('Barcodes','UMIs')],by='Barcodes')
      
      
      if (length(which(cross_comparison$UMIs > Test_BC_UMIs)) < 1) {
        BCs_to_consolidate_UMI <- rbind(BCs_to_consolidate_UMI,data_with_EdDistCnt[i,])
        percent_of_parent_umi <- cross_comparison
        colnames(percent_of_parent_umi) <- c('Barcodes','Pct_of_Parent')
        
        cross_comparisons[[data_with_EdDistCnt$Barcodes[i]]] <- cross_comparison
        for (x in 1:nrow(cross_comparison)) {
          percent_of_parent_umi$Pct_of_Parent[x] <- (cross_comparison$UMIs[x]/Test_BC_UMIs)*100
          cross_comparison$UMIs[x] <- Test_BC_UMIs - cross_comparison$UMIs[x]
        }
        avg_dist <- mean(cross_comparison$UMIs)
        #add_data <- data.frame(Barcodes = data_with_EdDistCnt[i,c('Barcodes')],Avg_UMI_Diff = as.numeric(avg_dist))
        #Test_BCcons_UMI_Dist <- rbind(Test_BCcons_UMI_Dist,add_data)
        
        BCcons_Indiv_Comps[[data_with_EdDistCnt$Barcodes[i]]] <- percent_of_parent_umi
        
      }
    }
  }
  
  
  if (length(BCcons_Indiv_Comps) > 0) {
    
    for (i in length(BCcons_Indiv_Comps):1) {
      length_check <- nrow(filter(BCcons_Indiv_Comps[[i]], Pct_of_Parent <= consolidation_percentage))
      if (length_check > 0) {
        BCcons_Indiv_Comps[[names(BCcons_Indiv_Comps)[i]]] <- filter(BCcons_Indiv_Comps[[i]], Pct_of_Parent <= consolidation_percentage)
      } else {
        BCcons_Indiv_Comps <- BCcons_Indiv_Comps[-i]
      }
    }
    
    
    ####################################################################################
    ##
    ##    Collapse edit distance barcodes into their parent barcode
    ##
    ####################################################################################
    
    ## I'm going to take advantage of that summarize function to do the barcode collapse.
    #     But instead of that weird thing the AI was doing, I'm going to make a clusters
    #     column and use it correctly.
    BCs_to_consolidate <- data.frame(Barcodes = names(BCcons_Indiv_Comps))
    BCs_to_consolidate <- left_join(BCs_to_consolidate,BCs_to_consolidate_UMI,by='Barcodes')
    BCs_to_consolidate <- BCs_to_consolidate[order(BCs_to_consolidate$UMIs),]
    
    
    ## I created this new copy of the data variable for testing purposes. It is unnecessary,
    #     but it's in the code now. It wouldn't be too hard to remove it, however.
    consolidated_data <- data
    
    ## I'm adding a column called "clusters" to the input dataset. At first, it is just
    #     assigning each barcode to it's own unique cluster
    consolidated_data$clusters <- c(1:nrow(data))
    
    ## This for loop here is going to identify all the barcodes within edit_distance
    #     of the 'parent' barcodes and assign them and their parent barcode new cluster values
    #  The first of these new clusters is assigned the value of nrow(data) + 1 so we can be
    #     sure we are using cluster values that haven't been used yet.
    cluster_ID <- nrow(data) + 1
    
    ## NOTE: If you wanted to use a different 'parent barcode' list for consolidation. You would
    #     replace the BCs_to_consolidate_UMI variable in the line of code below.
    for (i in BCs_to_consolidate$Barcodes) {
      
      ## Get the index of the parent barcode
      bcs_within_edit <- grep(i,data$Barcodes)
      
      for (x in 1:nrow(BCcons_Indiv_Comps[[i]])) {
        
        bcs_within_edit[length(bcs_within_edit)+1] <- grep(BCcons_Indiv_Comps[[i]]$Barcodes[x],data$Barcodes)
        
      }
      
      ## Get the indices of the other barcodes within edit_distance of the parent barcode
      #bcs_within_edit <- which(pairs[,2] == BC_ind)
      #bcs_within_edit <- pairs[bcs_within_edit,1]
      
      ## Add the parent barcode's index to this list so it gets assigned to the same
      #     cluster as it's associated barcodes
      #bcs_within_edit[length(bcs_within_edit)+1] <- BC_ind
      
      #print(paste0(as.character(length(which(consolidated_data$clusters[bcs_within_edit] > nrow(data))))))
      
      ## Assign a unique cluster ID to this list of barcodes
      consolidated_data$clusters[bcs_within_edit] <- cluster_ID
      
      ## Make sure the next cluster will have a unique cluster ID
      cluster_ID <- cluster_ID + 1
      
    }
    
    ## In case the data isn't ordered by UMI counts, we need to make sure it is because
    #     the summarise function is going to collapse the data into the entry that is "first"
    #     in the data being fed into it. So, we are making sure our highest UMI count barcode will be
    #     higher in the list than any of it's associated barcodes.
    consolidated_data <- consolidated_data[order(consolidated_data$UMIs,decreasing=TRUE),]
    
    ##  Now we call the summarise function telling it that the groups we want summarized are defined by
    #     the values in the cluster column we just created.
    consolidated_data <- consolidated_data %>%
      group_by(clusters) %>%
      summarise(Barcodes = first(Barcodes),
                Peptide = first(Peptide),
                RawCounts = sum(RawCounts, na.rm = TRUE),
                UMIs = sum(UMIs, na.rm = TRUE),
                pctTotalReads = sum(pctTotalReads, na.rm = TRUE))
    
    consolidated_data <- as.data.frame(consolidated_data[,c('Barcodes','Peptide','RawCounts','UMIs','pctTotalReads')])
    
    ##  These variables are just a reality check to ensure that the consolidated data was reduced by the
    #       number of barcodes that we would expect.
    Consolidated_BC_Total <- sum(BCs_to_consolidate_UMI$BCs_within_edit_dist)
    Predicted_Consolidated_DataFrameSize <- nrow(data) - Consolidated_BC_Total
    Actual_Consolidated_DataFrameSize <- nrow(consolidated_data)
    
    print(paste0('Original Data Size = ',as.character(nrow(data)),
                 ' rows. Pre Pct-of-Parent Filter Output Size = ',as.character(Predicted_Consolidated_DataFrameSize),
                 ' rows. Actual Output Size = ',as.character(Actual_Consolidated_DataFrameSize),' rows.'))
    
    ## Adding back in the UMI Rankings column and adding back the original column names
    consolidated_data <- consolidated_data[order(consolidated_data$UMIs, decreasing = TRUE),]
    colnames(consolidated_data) <- c('Barcodes','Peptide','RawCounts','UMICounts','pctTotalReads')
    
    #consolidated_data <- consolidated_data[order(consolidated_data$UMICounts,decreasing = TRUE),]
    
    return(consolidated_data)
    
  } else {
    
    print(paste0('Original Data Size = ',as.character(nrow(data)),
                 ' rows. No sequencing error barcodes identified in this dataset'))
    
    colnames(data) <- c('Barcodes','Peptide','RawCounts','UMICounts','pctTotalReads')
    return(data)
  }
  
}









