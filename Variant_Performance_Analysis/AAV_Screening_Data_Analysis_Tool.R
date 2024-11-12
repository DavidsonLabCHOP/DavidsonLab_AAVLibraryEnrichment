# This script is set up to process Illumina AAV Amplicon Screening Data based on
#     experiments performed by the Davidson Lab

# NOTE: The nucleic acid sequence encoding the random hepatamer peptide inserts in
#       the PM-AAVs is referred to throughout this script using the term 'Barcode'


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
library(pheatmap)

## This toggle will cause the script to delete all the variables it creates once it has finished
#     running and free up any unused memory.
## CAUTION: Setting this variable to TRUE and then running this script to completion will
#     delete all variables in your Environment (not just the ones created by this script).
Clear_Environment_Upon_Completion <- FALSE


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
UMI_Per_Barcode_Threshold <- 100




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
Save_File_Tag <-  "Your_Experiment_ID"

## This option will save *.rds objects containing the lists of UMI counts and Enrichment metric tables created by
#     this script.
Save_Results_Tables <-  TRUE

## Toggle determining whether the code will save heatmaps showing the top performing AAV variants across the set of 
#     experimental data being analyzed.  
#  Separate heatmaps will be saved for the UMI Counts metric and the Enrichment metric. These plots will show data for the
#     top <n_BCs_to_plot> barcodes (i.e., variants) as determined by average metric ranking across the experimental datasets.
#  NOTE1a: If too many barcodes are plotted (problems tend to begin around 30+), the peptide sequences may overlap and be difficult to read
#  NOTE1b: Similarly, if too many experimental datasets are present, their labels may also end up overlapping.
#  NOTE2: Each dataset will be labeled using its filename (sans "-counts.txt"). If your filenames are too long, they will
#             cause the resulting heatmap to be relatively small and some of the dataset names may get clipped and not
#             display fully.
Save_Region_Comparison_Heatmaps <- FALSE

## Toggle determining whether the code will save PDFs containing bar plots of the top hits for each dataset. Separate bar plots will be 
#     created for UMI Counts and Enrichment metrics.
## NOTE: A 'UMI Counts - Top Hits' and 'Enrichment - Top Hits' bar plot will be created for EACH experimental dataset. So, if you
#     have 30 datasets, this toggle will end up creating 60 PDF files.
Save_Dataset_TopHits_Barplots <- FALSE

## This variable determines how many barcodes (variants) are plotted in the analysis figures created if the 
#     Save_Region_Comparison_Heatmaps and/or Save_Dataset_TopHits_Barplots are set to TRUE. 
#     It is recommended to keep this variable at 30 or less to prevent overlap in peptide labeling.
n_BCs_to_plot <- 25


## The following code determines where the output of this script is stored
## Change Save_Folder_Parent_Dir and/or Save_Dir to control where the results tables are saved. 
#     Currently, the script creates and outputs to a folder called 'Results-<Save_File_Tag>' 
#     which it creates in Save_Folder_Parent_Dir
if (Save_Results_Tables | Save_Region_Comparison_Heatmaps | Save_Dataset_TopHits_Barplots) {
  Save_Folder_Parent_Dir <- getwd()
  Save_Dir = paste0(Save_Folder_Parent_Dir,"/Results","-",Save_File_Tag)
  dir.create(Save_Dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
}






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
AAV_Input_Vector_ID <- 'Name_of_your_input_vector_counts_txt_file' 







  
  





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
## Create a function to import datasets  
importdatasets <- function(sample) {
  outputDF <- data.frame(read_table(sample, col_names = FALSE ))
  return(outputDF)
}

### Apply the importdatasets function to import the datasets read from the text document
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




#########################
#################################################
##  Apply UMI_Per_Barcode Threshold filter

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
      bcs_within_edit <- grep(i,dataset$Barcodes)
      
      ## Get the index of each barcode in the <edit_dist_setting> group that was identified as
      #     a sequencing error due to failing the Pct_of_Parent_Threshold check.
      for (x in 1:nrow(BCs_within_edit_of_parent[[i]])) {
        bcs_within_edit[length(bcs_within_edit)+1] <- grep(BCs_within_edit_of_parent[[i]]$Barcodes[x],dataset$Barcodes)
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

# This code creates a reference for every barcode and peptide sequence observed across all datasets.
#     The reason this is necessary is because the same peptide sequence can be generated by multiple 
#     nucleic acid sequences (i.e., barcodes). This 'Final_Key' creates a unique name for each
#     peptide sequence (by appending #'s onto the duplicates). In addition to its analytical utility,
#     this key also serves a practical purpose, as many plotting functions require unique values
#     for categorical labels.
# This key will be saved at the end of the script if Save_Results_Tables is set to TRUE.

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
#     The second data frame contains Enrichment metric values, along with the percent_UMI
#     values for experimental and input vector data that are used in the Enrichment calculation.
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




###################################
##  Run Datasets                ###
###################################

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
  processed_data <- Create_Metric_Data_Frames(input_vector_dataset, exp_data)
  counter <<- counter + 1
  return(processed_data)
}

Metric_Tables <- lapply(alldatasets, RunAAV_Samples)

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





#########################################################################################
#########################################################################################
##
##    Output Save Options
##    
##    Note: Heatmaps and Bar Plots are only generated if Save_Region_Comparison_Heatmaps
##          and/or Save_Dataset_TopHits_Barplots is set to TRUE. You will find the 
##          functions to generate those plots under the relevant toggles below.
##
#########################################################################


setwd(Save_Dir)

# This variable is used in the Settings_And_Info output file
saved_files <- array()
saved_files[1] <- "Names of Saved Datafiles:"
saved_files[2] <- "----------------------"


## Save the enrichment and UMI counts data frames created by RunAAV_Samples
if (Save_Results_Tables) {
  
  saveRDS(Enrichment_Tables, file = paste0(as.character(Save_Dir),"/Enrichment_Tables-",Save_File_Tag,".rds"))
  saveRDS(UMIcount_Tables, file = paste0(as.character(Save_Dir),"/UMIcount_Tables-",Save_File_Tag,".rds"))
  saveRDS(AAV_Final_Key, file = paste0(as.character(Save_Dir),"/AAV_Final_Key-",Save_File_Tag,".rds"))
  
  saved_files[length(saved_files)+1] <- paste0("Enrichment_Tables-",Save_File_Tag,".rds")
  saved_files[length(saved_files)+1] <- paste0("UMIcount_Tables-",Save_File_Tag,".rds")
  saved_files[length(saved_files)+1] <- paste0("AAV_Final_Key-",Save_File_Tag,".rds")
  
}

## Save region comparison heatmaps and/or top hits bar plots
if (Save_Region_Comparison_Heatmaps | Save_Dataset_TopHits_Barplots) {
  
  ##############################################################################
  ## This function combines metric data for the various experimental datasets
  #     into a single data frame to be used for plotting
  combine_datasets <- function(data_list, metric) {
    
    ## Removing the input vector from the list of experimental datasets
    input_index <- grep(AAV_Input_Vector_ID, names(data_list))
    data_list <- data_list[-input_index]
    
    ## Creating a master list of all barcodes present in experimental datasets.
    #     Then cross-referencing it against the AAV_Final_Key data frame as a single
    #     source of name resolution for duplicate peptide sequences.
    combined_data <- do.call(rbind, data_list)
    combined_data <- combined_data[!duplicated(combined_data$Barcodes),]
    combined_data <- data.frame(Barcodes = combined_data$Barcodes)
    combined_data <- left_join(combined_data, AAV_Final_Key, by='Barcodes')
    
    if (metric == 'Enrichment') {
      for (i in 1:length(data_list)) {
        combined_data <- left_join(combined_data, data_list[[i]][,c('Barcodes','Enrichment')], by='Barcodes')
        colnames(combined_data)[ncol(combined_data)] <- gsub('-counts.txt','',names(data_list)[i])
      }
      combined_data[is.na(combined_data)] <- 0
      combined_data$AVG <- rowMeans(combined_data[,c(3:ncol(combined_data))])
      combined_data <- relocate(combined_data, AVG, .after = Peptide)
      combined_data <- arrange(combined_data, desc(AVG))
      return(combined_data)
      
    } else if (metric == 'UMICounts') {
      for (i in 1:length(data_list)) {
        combined_data <- left_join(combined_data, data_list[[i]][,c('Barcodes','UMI_Exp')], by='Barcodes')
        colnames(combined_data)[ncol(combined_data)] <- gsub('-counts.txt','',names(data_list)[i])
      }
      combined_data[is.na(combined_data)] <- 0
      combined_data$AVG <- rowMeans(combined_data[,c(3:ncol(combined_data))])
      combined_data <- relocate(combined_data, AVG, .after = Peptide)
      combined_data <- arrange(combined_data, desc(AVG))
      return(combined_data)
      
    } else {
      print('Invalid metric argument')
      return(NULL)
    }
  }
  
  ## Heatmap plotting function
  plot_region_comparison_heatmap <- function(Heat_input, n_BCs=30) {
    if (nrow(Heat_input>20)) {
      rowfontsize = 8
    } else {
      rowfontsize = 10
    }
    Heat_input[Heat_input==0] <- NA
    find_data <- grep('AVG',colnames(Heat_input))
    Heat_input <- dplyr::arrange(Heat_input, desc(AVG))
    if (nrow(Heat_input) < n_BCs) {
      n_BCs <- nrow(Heat_input)
    }
    Heat_input <- Heat_input[c(1:n_BCs),]
    rownames(Heat_input) <- Heat_input$Peptide
    
    output_graph <- pheatmap(Heat_input[,c((find_data+1):ncol(Heat_input))],        
                             labels_row = Heat_input$Peptide,
                             cluster_cols = FALSE,
                             #color = c.pal,
                             cluster_rows = FALSE,
                             dendrogram="none",
                             trace="none",
                             scale = "column",
                             fontsize = rowfontsize,
                             margins = c(2, 2),
                             border_color=NA,
                             na_col = "white",
                             angle_col = 315,
                             #gaps_row = 27,
                             fontsize_row = rowfontsize)
    return(output_graph)  
  }
  
  ## Tophits barplot plotting function
  plot_dataset_tophits_barplot <- function(input_dataset, metric_label, n_BCs=30) {
    barplot_list <- list()
    
    ## Creating a barplot for each dataset present in the input_dataset
    for (i in 4:ncol(input_dataset)) {
      bar_data <- input_dataset[,c(1,i)]
      colnames(bar_data) <- c('Barcodes','Metric')
      bar_data <- dplyr::arrange(bar_data, desc(Metric))
      
      ## Pulling Peptide labels from the AAV_Final_Key to ensure consistent names are used
      #     for any duplicate peptide sequences.
      bar_data <- left_join(bar_data, AAV_Final_Key, by='Barcodes')
      if (nrow(bar_data) < n_BCs) {
        n_BCs <- nrow(bar_data)
      }
      bar_data <- bar_data[c(1:n_BCs),]
      bar_data <- dplyr::arrange(bar_data, Metric)
      bar_data$Peptide <- factor(bar_data$Peptide, levels = bar_data$Peptide)
      
      Bar_plot_output <- ggplot(bar_data, aes(x=Peptide, y=Metric)) + 
        geom_col(aes(Metric,Peptide), color = "black", fill = "skyblue4", linewidth = 0.75) + 
        theme_classic() +
        scale_x_continuous(expand = c(0,0), position = "top") +
        labs(title = as.character(colnames(input_dataset)[i])) +
        theme(plot.title = element_text(hjust = 0.5, vjust = 3)) +
        theme(axis.title.y=element_blank()) +
        labs(x = as.character(metric_label))
      
      barplot_list[[colnames(input_dataset)[i]]] <- Bar_plot_output
    }
    return(barplot_list)
  }
  
  ## Function to save heatmap as a PDF
  save_heatmap_pdf <- function(x, filename, width=7, height=7) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  ## Function to save bar plot as a PDF
  save_barplot_pdf <- function(x, filename, width=7, height=7) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    invisible(print(x))
    dev.off()
  }
  
  
  ## Create data frames that will be used as input to the plotting functions
  Enrichment_dataframe <- combine_datasets(Enrichment_Tables, 'Enrichment')
  UMIcount_dataframe <- combine_datasets(UMIcount_Tables, 'UMICounts')
  
  
  ## Create experimental dataset (region) comparison heatmaps
  if (Save_Region_Comparison_Heatmaps) {
    
    ## Get plots from heatmap plotting functions
    region_comparison_heatmap_plot_Enrichment <- plot_region_comparison_heatmap(Enrichment_dataframe, n_BCs_to_plot)
    region_comparison_heatmap_plot_UMIs <- plot_region_comparison_heatmap(UMIcount_dataframe, n_BCs_to_plot)
    
    ## Save plots as PDFs
    save_filename <- paste0(Save_Dir,'/RegionComp-Enrichment-Heatmap-',Save_File_Tag,'.pdf')
    save_heatmap_pdf(region_comparison_heatmap_plot_Enrichment, save_filename)
    saved_files[length(saved_files)+1] <- save_filename
    
    save_filename <- paste0(Save_Dir,'/RegionComp-UMIs-Heatmap-',Save_File_Tag,'.pdf')
    save_heatmap_pdf(region_comparison_heatmap_plot_UMIs, save_filename)
    saved_files[length(saved_files)+1] <- save_filename
  }
  
  ## Create top hit bar plot for each experimental dataset
  if (Save_Dataset_TopHits_Barplots) {
    
    ## The barplot function outputs a list of plots, one for each dataset
    tophit_barplots_list_Enrichment <- plot_dataset_tophits_barplot(Enrichment_dataframe, 'Enrichment', n_BCs_to_plot)
    tophit_barplots_list_UMIs <- plot_dataset_tophits_barplot(UMIcount_dataframe, 'UMI Counts', n_BCs_to_plot)
    
    ## Output the lists of bar plots as PDFs
    for (i in 1:length(tophit_barplots_list_Enrichment)) {
      save_filename <- paste0(Save_Dir,'/TopHits-Enrichment-',names(tophit_barplots_list_Enrichment)[i],'.pdf')
      save_barplot_pdf(tophit_barplots_list_Enrichment[[i]], save_filename)
      saved_files[length(saved_files)+1] <- save_filename
    }
    for (i in 1:length(tophit_barplots_list_UMIs)) {
      save_filename <- paste0(Save_Dir,'/TopHits-UMIs-',names(tophit_barplots_list_UMIs)[i],'.pdf')
      save_barplot_pdf(tophit_barplots_list_UMIs[[i]], save_filename)
      saved_files[length(saved_files)+1] <- save_filename
    }
  }
}





################################################################################
##  Generate Seetings_And_Info text file containing run information          ###
################################################################################

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
print("Processing and Filtering Parameters")
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
print(paste("Save_Region_Comparison_Heatmaps =", Save_Region_Comparison_Heatmaps))
print(paste("Save_Dataset_TopHits_Barplots =", Save_Dataset_TopHits_Barplots))
if (Save_Region_Comparison_Heatmaps | Save_Dataset_TopHits_Barplots) {
  print(paste("Number_of_BCs_to_Plot =", n_BCs_to_plot))
}
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



## CAUTION: Running this line of code (assuming Clear_Environment_Upon_Completion = TRUE) will
#     delete all variables in your Environment (not just the ones created by this script).
if (Clear_Environment_Upon_Completion) {rm(list=ls()); gc()}

