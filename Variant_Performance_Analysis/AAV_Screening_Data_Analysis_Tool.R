# This script is set up to process Illumina AAV Amplicon Screening Data based on
#     experiments performed by the Davidson Lab

# NOTE on Nomenclature: The nucleic acid sequence encoding the random heptamer peptide 
#     inserts in the PM-AAVs is referred to throughout this script using the term "Barcode"

################################
##  Loading Dependencies      ##
################################
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


## This toggle will cause the script to delete all the variables it creates once it has finished
#     running and free up any unused memory. 
Clear_Environment_Upon_Completion <- TRUE


##############################################
##  Important code settings and variables   ##
##############################################

#####################################################
##  Path to Data File List                        ###
#####################################################

## The only input required to run this script is the full path to a text file (UFT-8 encoded) containing the list
#     of data files (counts.txt output from the AAV_Barcode_Parsing_Tool.py) to be analyzed. Each data file name
#     in a separate row in the text file with no headers.
#     This text file must be located in the same directory as the experimental data. The script will set the
#     working directory to dirname(Data_File); so, the text file only needs to contain the file names, not full paths.
Data_File <- "Full_Path_to_Text_File_Containing_List_of_Counts_Files_to_be_Analyzed"

#Data_File <- "C:/Users/lewandowsb/OneDrive - Children's Hospital of Philadelphia/Desktop/Broad ICV Analysis/Round 2/Testset_Temp/SampleIDs.txt"

Data_File_List <- basename(Data_File)
Experimental_Data_Dir <-  dirname(Data_File)
setwd(Experimental_Data_Dir)



#######################################
###  UMI Counts Threshold           ###
#######################################

## This variable defines the UMI count threshold used to filter out weakly detected Barcodes.
#     This value indicates the minimum number of UMIs required for a given Barcode to be included 
#     in enrichment analyses. This threshold will need to be determined empirically based on the
#     distribution of values across a given set of screening datasets.
UMI_Per_Barcode_Threshold <- 10




##########################################################
##  Sequencing Error Detection and Removal Options     ###
##########################################################

## This toggle determines whether the script will perform automated sequence error barcode detection
##    and removal. Note: This is a memory intensive process. It is run individually on each dataset. 
##    If a given dataset includes too many barcodes, it is possible that R could run out of available memory. 
##    The UMI_Per_Barcode_Threshold is applied before the sequence error detection/removal is run. So, you can
##    try increasing that threshold if you are having memory issues.
Collapse_Barcodes <- TRUE

## This variable determines the edit distance between barcodes used to determine which barcodes
##    are compared as potential sequencing errors.
edit_dist_setting <- 1

## Putative sequencing error barcodes need to have less than this Percent_of_Parent_Threshold
#     in order to be considered an actual sequencing error barcode. Percent of Parent is calculated
#     using the following equation for each putative sequencing error barcode
#     (Err_Barcode_UMIs / Parent_Barcode_UMIs) * 100
#     Note: this value is a % (i.e., it ranges from 0 to 100, not 0 to 1).
Percent_of_Parent_Threshold <- 1





###############################
###  Script Output Options  ###
###############################

## The string saved to the Save_File_Tag variable will be appended to the end of the standard names for each of the outputs
#     saved by this script. You'll want to make this unique to each data set so you don't accidentally over-write
#     the output of a previous run.
Save_File_Tag = "Your_Experiment_ID"

## This option will save *.rds objects containing the lists of UMI counts and Enrichment metric tables created by
#     this script.
Save_Results_Tables = TRUE

## Change Save_Parent_Dir and/or Save_Dir to control where the results tables are saved. 
#     Use getwd() to save in same location as the experimental files.
#     Currently, the script creates and outputs to a folder called 'Results-<Save_File_Tag>' in Save_Parent_Dir 
Save_Parent_Dir <- getwd()
Save_Dir = paste0(Save_Parent_Dir,"/Results","-",Save_File_Tag)
dir.create(Save_Dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")





#########################################
###  Identification of Input Vectors  ###
#########################################

### Normalization and data analysis in this script requires that the input vector files be analyzed separately 
#     from the experimental results files. You can enter the Input Vector file IDs manually or try an automated
#     identification of input files by setting AutoFind_Input_Files = TRUE.
#     The auto-find feature relies on the Input Vector files having "Input" in their file name.
AutoFind_Input_Files <- TRUE

# If your input vector's data file name does not contain Input, you will need to set AutoFind_Input_File to FALSE
#    and then enter the file name of your input vector below. This file name must match the entry for the input vector
#    data file in the list identified by the Data_File variable (defined above)
# If AutoFind_Input_Files <- TRUE, then it doesn't matter what is entered here; it will be over-written later in the script.
AAV_Input_Vector_ID <- 'Name of your input vector counts.txt file' 







  
  





##################################
###############################################################
################################################################################
##
##  Begin Automated Screening Data Analysis
##
################################################################################
##
##  If the variables above are properly set, the code below should run without
##      need for further user input.
##
################################################################################
###############################################################
##################################



###############################################################
###   Setting up Run Information Tracking Variables         ###
###############################################################

# Setting up variables that will track certain activity in the script
#    (e.g., removal of problematic datasets).
# The information they store is output as part of the Settings_and_Info text 
#    file that is saved at the end of the script.
QC_Report <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(QC_Report) <- c("Data_File","Reason_For_Exclusion")

Output_Notes <<- array()
Output_Notes[1] <- "Additional Information:"
Output_Notes[2] <- "-----------------------"





####################################
#########################################################
##
## Import the counts data files created by running AAV_Barcode_Parsing_Tool.py
##   on amplicon FASTQ files

# These *counts.txt files contain four columns of data:
# Column 1 is the Barcode sequence in DNA base pairs
# Column 2 is the Barcode sequence translated to amino acid letter code
# Column 3 is the total counts
# Column 4 is the UMI normalized counts


### Create a list of sample IDs
SampleIDs <- data.frame(read_csv(Data_File_List, col_names=FALSE))
names(SampleIDs) <- c("Names") 



##################################################################
## Create a function to import the dataset of interest        

importdatasets <- function(sample) {
  outputDF <- data.frame(read_table(sample, col_names = FALSE ))
  return(outputDF)
}

### Apply the importdatasets function to read in the datasets read from the text document
#       indicated by the Data_File variable
alldatasets=list()
alldatasets <- lapply(SampleIDs$Names, FUN = importdatasets)
names(alldatasets) <- SampleIDs$Names  

for (i in c(1:length(alldatasets))) {
  colnames(alldatasets[[i]]) <- c("Barcodes","Peptide","RawCounts","UMICounts")
}











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

# Adding input vector file ID info to the run tracking variable
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



#########################################################################################
#########################################################################################
##
##    COLLAPSE BARCODE FUNCTION
##
##    The toggle to control whether this function is run can be found at the beginning
##        of the script. In addition, the variables that determine the error rate (edit_dist_setting)
##        and the Pct_of_Parent_Threshold can also be found at the top of the script.
##
##    This function will find barcodes within edit_dist_setting of each other and combine
##        then into the a single entry based on which one of them has the highest UMI counts.
##
#########################################################################


collapse_BCs <- function(dataset) {
  
  ###############################################
  ##  Calculate edit_distances and identify oPool barcode pairs
  ##     within <edit_dist_setting> of each other.
  
  ## Calculate the pairwise edit distance between strings
  #  If you are looking for the step that takes the longest to run in this code, this is it. The adist()
  #    function creates a pair-wise matrix which grows exponentially with the number of comparisons (Barcodes).
  #  The adist() has an argument (costs) to define the costs of various string manipulations.
  #    We are only interested in substitutions. To prevent insertion or deletion based results registering as
  #    valid barcode comparisons, the cost of insertions and deletions is set to 10.
  #        In contrast, the cost of substitutions is set to 1, so the value of the edit_dist_setting variable  
  #        will be directly correlated with the maximum allowable number of edits.
  edit_costs <- list()
  edit_costs[['insertions']] <- 10
  edit_costs[['deletions']] <- 10
  edit_costs[['substitutions']] <- 1
  
  edit_distances <- adist(dataset$Barcodes,costs=edit_costs,counts=TRUE)
  
  # Filter the pairs with an edit distance of <edit_dist_setting> or less (excluding the diagonal elements)
  pairs <- which(edit_distances <= edit_dist_setting & edit_distances > 0, arr.ind = TRUE)
  
  
  #################################################
  ##  Creating a new data frame, barcode_similarities, to store the count of other
  #      barcodes within <edit_dist_setting> of each individual barcode. This data frame is a simple 
  #      list of the barcodes in 'dataset' and a column called 'BCs_within_edit_dist'.
  #   It is not strictly necessary to create a separate data frame for this purpose, but
  #      the barcode_similarities data frame can be a useful troubleshooting tool if you find
  #      your data behaving strangely.
  
  ## Create a new empty data frame and populate it with the barcode info from 'dataset'
  barcode_similarities <- data.frame(matrix(nrow=nrow(dataset),ncol=2))
  colnames(barcode_similarities) <- c('Barcodes','BCs_within_edit_dist')
  barcode_similarities$Barcodes <- dataset$Barcodes
  
  for (i in 1:nrow(dataset)) {
    ## Count barcodes within <edit_distance_setting> of barcode_similarities$Barcodes[i]
    barcode_similarities$BCs_within_edit_dist[i] <- length(which(pairs[,2] == i))
  }
  
  ## Adding the BCs_within_edit_dist information to the dataset variable
  dataset <- cbind(dataset,barcode_similarities$BCs_within_edit_dist)
  colnames(dataset) <- c('Barcodes','Peptide','RawCounts','UMICounts','BCs_within_edit_dist')
  rm(barcode_similarities)
  
  #####################################################################################
  ##
  ##  Create a data frame containing the "parent" barcodes that should be used as the 
  ##      primary input that the other barcodes will be collapsed into
  ##
  ##  Parent barcodes are defined as the barcode with the greatest UMI counts amongst a
  ##      group of barcodes within <edit_dist_setting> of each other
  ##
  ######################################################################################
  
  ## Creating a list structure to store sets of barcodes within <edit_dist_setting> of each other
  #     that might contain sequencing errors
  BCs_within_edit_of_parent <- list()
  
  ## Creating an empty data frame to store the identity of barcodes identified as Parent barcodes
  #     by the subsequent for loop
  Putative_Parent_BCs <- data.frame(matrix(nrow=0,ncol=5))
  colnames(Putative_Parent_BCs) <- c('Barcodes','Peptide','RawCounts','UMICounts','BCs_within_edit_dist')
  
  ## This loop will search through the data and identify the parent barcodes for each set of
  #     barcodes within <edit_dist_setting> of each other. It will save this information to the
  #     variables created above. In the following step, this information will be used to
  #     identify parent barcode / sequencing error pairs.
  
  # For each barcode in the current dataset
  for (i in 1:nrow(dataset)) {
    
    # Check if current barcode is within <edit_dist_setting> of any other barcodes
    if (dataset$BCs_within_edit_dist[i] > 0) {
      
      # UMI counts of the barcode currently being tested
      Test_BC_UMIs <- dataset$UMICounts[i]
      
      # Grabbing indices of other barcodes within 1 edit distance of the barcode currently being tested
      #    NOTE: It is critical that dataset is not re-ordered between calculating the pairs matrix
      #             and running this code. Otherwise, the indices identified here will not map to the
      #             correct barcodes in dataset.
      bcs_within_edit <- which(pairs[,2] == i)
      bcs_within_edit <- pairs[bcs_within_edit,1]
      
      ## Pulling information for the current set of barcodes being analyzed for putative sequencing errors
      cross_comparison <- dataset[bcs_within_edit, c('Barcodes','UMICounts')]
      
      ## Check if the current barcode is the "parent barcode" within the current set of putative sequencing errors.
      #     The parent barcode is defined as the barcode within the set with the highest UMI count.
      if (length(which(cross_comparison$UMICounts > Test_BC_UMIs)) < 1) {
        
        ## If barcode has the highest UMI count amongst its set, store it in the putative parent barcode data frame.
        Putative_Parent_BCs <- rbind(Putative_Parent_BCs,dataset[i,])
        
        ## For all the other barcodes in the <edit_dist_setting> group, a % of parent UMIs metric is calculated
        #     Pct_of_Parent = ([Barcode UMI count] / Parent Barcode UMI count) * 100
        cross_comparison$Pct_of_Parent <- (cross_comparison$UMICounts/Test_BC_UMIs)*100
        cross_comparison$UMICounts <- NULL
        
        BCs_within_edit_of_parent[[dataset$Barcodes[i]]] <- cross_comparison
        
      }
    }
  }
  
  
  ##################################################################
  ##
  ##  Test for Parent Barcode / Sequencing Error pairs
  
  ## Test if any parent barcode sequences were detected in the current dataset
  if (length(BCs_within_edit_of_parent) > 0) {
    
    ## For each element in the BCs_within_edit_of_parent list working backwards
    for (i in length(BCs_within_edit_of_parent):1) {
      ## Test if any barcodes in the <edit_dist_setting> list have Pct_of_Parent values below the <Percent_of_Parent_Threshold> threshold.
      #     Barcodes with Pct_of_Parent <= <Percent_of_Parent_Threshold> are considered sequencing errors
      length_check <- nrow(dplyr::filter(BCs_within_edit_of_parent[[i]], Pct_of_Parent <= Percent_of_Parent_Threshold))
      if (length_check > 0) {
        ## Filter the list of <edit_dist_setting> barcodes to remove any that passed the Pct_of_Parent threshold check, leaving only the
        #     barcodes that failed the check and are considered sequencing errors.
        BCs_within_edit_of_parent[[names(BCs_within_edit_of_parent)[i]]] <- dplyr::filter(BCs_within_edit_of_parent[[i]], Pct_of_Parent <= Percent_of_Parent_Threshold)
      } else {
        ## If no barcodes within the <edit_dist_setting> group fail the Pct_of_Parent threshold check (i.e., none are putative sequencing errors),
        #     then remove the Parent Barcode entry from the BCs_within_edit_of_parent list. We don't need to retain Parent Barcode entries
        #     that are associated with no putative sequencing errors.
        BCs_within_edit_of_parent <- BCs_within_edit_of_parent[-i]
      }
    }
    
    
    ####################################################################################
    ##
    ##    Collapse putative sequencing error barcodes into their parent barcode
    ##
    ####################################################################################
    
    BCs_to_consolidate <- Putative_Parent_BCs[Putative_Parent_BCs$Barcodes %in% names(BCs_within_edit_of_parent),]
    BCs_to_consolidate <- dplyr::arrange(BCs_to_consolidate, UMICounts)
    
    #BCs_to_consolidate <- data.frame(Barcodes = names(BCs_within_edit_of_parent))
    #BCs_to_consolidate <- left_join(BCs_to_consolidate,Putative_Parent_BCs,by='Barcodes')
    #BCs_to_consolidate <- BCs_to_consolidate[order(BCs_to_consolidate$UMICounts),]
    
    ## Adding a column called "clusters" to dataset to create groupings that can be fed
    #     into a dplyr summarize function. To begin, all barcodes will be assigned to their
    #     own unique cluster.
    dataset$clusters <- c(1:nrow(dataset))
    
    ## Next, the indices for each set of parent barcode / sequencing error barcodes are identified
    #     and assigned their own unique cluster value (cluster_ID). The highest cluster value in
    #     the dataset at this point = nrow(dataset). So, we'll start with the value:
    #     nrow(data) + 1 to begin labeling our parent / seq error clusters.
    cluster_ID <- nrow(data) + 1
    
    ## Running through each parent barcode in our BCs_to_consolidate data frame
    for (i in BCs_to_consolidate$Barcodes) {
      ## Get the index of the parent barcode
      bcs_within_edit <- grep(i,data$Barcodes)
      
      ## Get the index of each barcode in the <edit_dist_setting> group that was identified as
      #     a sequencing error due to failing the Pct_of_Parent_Threshold check.
      for (x in 1:nrow(BCs_within_edit_of_parent[[i]])) {
        bcs_within_edit[length(bcs_within_edit)+1] <- grep(BCs_within_edit_of_parent[[i]]$Barcodes[x],data$Barcodes)
      }
      
      ## Assign a unique cluster ID to this list of barcodes
      dataset$clusters[bcs_within_edit] <- cluster_ID
      
      ## Make sure the next cluster will have a unique cluster ID
      cluster_ID <- cluster_ID + 1
    }
    
    ## Pulling total BCs pre-collapse for output to the command line to help track processing progress
    nBCs_pre_collapse <- nrow(dataset)
    
    ## It is important to ensure the dataset is ordered by UMICounts in descending order.
    #     The summarise function is going to collapse the data into the entry that is "first"
    #     in the data being fed into it. So, make sure the highest UMI count barcode is
    #     higher in the list than any of it's associated barcodes.
    dataset <- dplyr::arrange(dataset, desc(UMICounts))
    
    ## Call the summarise function using the Parent BC / Seq Error BC groups (clusters) that
    #     were assigned above.
    dataset <- dataset %>%
      group_by(clusters) %>%
      summarise(Barcodes = first(Barcodes),
                Peptide = first(Peptide),
                RawCounts = sum(RawCounts, na.rm = TRUE),
                UMICounts = sum(UMICounts, na.rm = TRUE))
    
    dataset <- as.data.frame(dataset[,c('Barcodes','Peptide','RawCounts','UMICounts')])
    dataset <- dplyr::arrange(dataset, desc(UMICounts))
    
    print(paste0('Original Data Size = ',as.character(nBCs_pre_collapse),
                 ' barcodes. Data Size After SeqErr Removal = ',as.character(nrow(dataset)),' barcodes.'))
    
    return(dataset)
    
  } else {
    
    ## If no parent barcode / sequencing errors pairs are identified in the current dataset,
    ##    reformat to remove changes and return dataset
    print(paste0('Original Data Size = ',as.character(nrow(dataset)),
                 ' barcodes. No sequencing error barcodes identified in this dataset'))
    
    dataset <- dataset[,c('Barcodes','Peptide','RawCounts','UMICounts')]
    dataset <- dplyr::arrange(dataset, desc(UMICounts))
    
    return(dataset)
  }
  
}



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
#   Create Metric Data Tables
#
#   This function creates two sets of data frames for each dataset. One data frame contains
#     UMICounts values for the experimental dataset and the input vector (for reference).
#     The second data frame contains results of Enrichment metric calculations, and the
#     UMIpct values for experimental and input vector data that are used in the Enrichment calculation.
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
  
  ## Creating Enrichment metric data frame
  experimental_data_enrichment <- experimental_data[,c('Barcodes','Peptide')]
  ## Calculating Percent UMI values for Experimental and Input Vector data
  experimental_data_enrichment$PctUMI_Input <- experimental_data$InputVector / sum(experimental_data$InputVector)
  experimental_data_enrichment$PctUMI_Exp <- experimental_data$ExpData / sum(experimental_data$ExpData)
  ## Calculating Enrichment values
  experimental_data_enrichment$Enrichment <- (experimental_data_enrichment$PctUMI_Exp - experimental_data_enrichment$PctUMI_Input)/experimental_data_enrichment$PctUMI_Input
  experimental_data_enrichment$Enrichment <- experimental_data_enrichment$Enrichment + 1
  
  experimental_data_enrichment <- dplyr::arrange(experimental_data_enrichment, desc(Enrichment))
  experimental_data_enrichment$Enrichment_Rank <- c(1:nrow(experimental_data_enrichment))
  
  ## Creating UMI Counts metric data frame
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
##    Create and Save AAV Metric Tables
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



##  This counter is used to provide some minor progress update functionality
counter <- 1

##  Simple wrapper function to run the Create_Metric_Data_Frames function with a
#     command line output indicating the dataset being processed

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
print(paste("Percent_of_Parent_Threshold =",as.character(Percent_of_Parent_Threshold)))
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



if (Clear_Environment_Upon_Completion) {rm(list=ls()); gc()}
