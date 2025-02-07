#!/usr/bin/env python

import sys
import argparse

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-D', '--inputDir', required=True, help='Directory containing input files')
parser.add_argument('-F', '--inputFastq', required=True, help='Input Fastq File')
args = parser.parse_args()

###############################################################
##
##  Critical User-Defined Parameters
##
##  You need to provide values into the variables below that define the regions containing
##      flanking sequences, the peptide insert, and the UMI based on the structure
##      of your amplicon. These values will probably change depending on the
##      serotype you are analyzing, so you may want to create separate versions of this
##      code for processing different serotypes.


# The code will look for perfect string matches to two flanking sequences as a quality control
#       check to ensure the amplicon is properly formed.
# Use the variables below to enter the indices and nucleic acid sequence defining your flanking sequences 
flank_1_start_index = 0
flank_1_end_index = 5
flank_1_sequence = 'ATCGA'

flank_2_start_index = 100
flank_2_end_index = 105
flank_2_sequence = 'ATCGA'

# This script assuming a heptamer peptide insert, so you only need to provide the first index
#       of the insert's location in the amplicon. If your insert is larger or smaller than 7
#       amino acids, you will need to adjust these variables and the code that controls the output
#       of the script.
peptide_modification_start_index = 50
peptide_modification_end_index = peptide_modification_start_index + 21

# Provide the start and end indices defining the location of the unique molecular identifier (UMI)
#       in your amplicon.
UMI_start_index = 120
UMI_end_index = 130

###############################################################
##################################

# Make empty dictionaries
Reads_Dict = {}
Barcodes_Dict = {}

# Create Codon dictionary
Codon_Dict = {}
Codon_Dict["TTT"] = "F"
Codon_Dict["TTC"] = "F"
Codon_Dict["TTA"] = "L"
Codon_Dict["TTG"] = "L"
Codon_Dict["CTT"] = "L"
Codon_Dict["CTC"] = "L"
Codon_Dict["CTA"] = "L"
Codon_Dict["CTG"] = "L"
Codon_Dict["ATT"] = "I"
Codon_Dict["ATC"] = "I"
Codon_Dict["ATA"] = "I"
Codon_Dict["ATG"] = "M"
Codon_Dict["GTT"] = "V"
Codon_Dict["GTC"] = "V"
Codon_Dict["GTA"] = "V"
Codon_Dict["GTG"] = "V"
Codon_Dict["TCT"] = "S"
Codon_Dict["TCC"] = "S"
Codon_Dict["TCA"] = "S"
Codon_Dict["TCG"] = "S"
Codon_Dict["CCT"] = "P"
Codon_Dict["CCC"] = "P"
Codon_Dict["CCA"] = "P"
Codon_Dict["CCG"] = "P"
Codon_Dict["ACT"] = "T"
Codon_Dict["ACC"] = "T"
Codon_Dict["ACA"] = "T"
Codon_Dict["ACG"] = "T"
Codon_Dict["GCT"] = "A"
Codon_Dict["GCC"] = "A"
Codon_Dict["GCA"] = "A"
Codon_Dict["GCG"] = "A"
Codon_Dict["TAT"] = "Y"
Codon_Dict["TAC"] = "Y"
Codon_Dict["TAA"] = "STOP"
Codon_Dict["TAG"] = "STOP"
Codon_Dict["CAT"] = "H"
Codon_Dict["CAC"] = "H"
Codon_Dict["CAA"] = "Q"
Codon_Dict["CAG"] = "Q"
Codon_Dict["AAT"] = "N"
Codon_Dict["AAC"] = "N"
Codon_Dict["AAA"] = "K"
Codon_Dict["AAG"] = "K"
Codon_Dict["GAT"] = "D"
Codon_Dict["GAC"] = "D"
Codon_Dict["GAA"] = "E"
Codon_Dict["GAG"] = "E"
Codon_Dict["TGT"] = "C"
Codon_Dict["TGC"] = "C"
Codon_Dict["TGA"] = "STOP"
Codon_Dict["TGG"] = "W"
Codon_Dict["CGT"] = "R"
Codon_Dict["CGC"] = "R"
Codon_Dict["CGA"] = "R"
Codon_Dict["CGG"] = "R"
Codon_Dict["AGT"] = "S"
Codon_Dict["AGC"] = "S"
Codon_Dict["AGA"] = "R"
Codon_Dict["AGG"] = "R"
Codon_Dict["GGT"] = "G"
Codon_Dict["GGC"] = "G"
Codon_Dict["GGA"] = "G"
Codon_Dict["GGG"] = "G"

# Read in input .fastq and store lines as a python dictionary
with open(args.inputFastq, "r") as infile:
    line_ct = 0
    for line in infile:
        if (line_ct % 4 == 1):
            Reads_Dict[line_ct]=str(line[0:])
        line_ct += 1

# Process every read stored in the Reads_Dict
for read in Reads_Dict.values():
    # Check for perfect string match to flanking sequence #1
    if (read[int(flank_1_start_index):int(flank_1_end_index)] == str(flank_1_sequence)):
        # Check for perfect string match to flanking sequence #2
        if (read[int(flank_2_start_index):int(flank_2_end_index)] == str(flank_2_sequence)):
            # If read passes flanking sequence checks, extract the peptide variant nucleic acid sequence (aka, Barcode)
            BARCODE=str(read[int(peptide_modification_start_index):int(peptide_modification_end_index)])
            # Extract UMI
            UMI=str(read[int(UMI_start_index):int(UMI_end_index)])
            
            # Ensure there are no unassigned base calls in the BARCODE sequence
            if "N" not in str(BARCODE):
                
                # Check if BARCODE has been previously observed.
                if (str(BARCODE) not in Barcodes_Dict.keys()):
                    # If barcode has not been previously observed, create a dictionary entry using the
                    #    nucleic acid sequence as the key
                    Barcodes_Dict[str(BARCODE)]=[]
                    # Add the UMI to the BARCODE's key entry
                    Barcodes_Dict[str(BARCODE)].append(str(UMI))
                else:
                    # Add the UMI to the BARCODE's key entry
                    Barcodes_Dict[str(BARCODE)].append(str(UMI))

# This script is designed to be run using a bash wrapper script that saves the output of this function to a text file.
# The code below outputs the contents of the Barcodes_Dict in lines containing four pieces of information each separated by a space:
#    1) The nucleic acid sequence of the peptide insert (BARCODE)
#    2) The amino acid sequence of the peptide insert
#    3) The raw read count
#    4) The UMI collapsed read count
for key in Barcodes_Dict:
    AA1 = key[0:3]
    AA2 = key[3:6]
    AA3 = key[6:9]
    AA4 = key[9:12]
    AA5 = key[12:15]
    AA6 = key[15:18]
    AA7 = key[18:21]
    print( key[:36] + " " + Codon_Dict[AA1] + Codon_Dict[AA2] + Codon_Dict[AA3] + Codon_Dict[AA4] + Codon_Dict[AA5] + Codon_Dict[AA6] + Codon_Dict[AA7] + " " + str(len(Barcodes_Dict[key])) + " " + str(len(list( dict.fromkeys(Barcodes_Dict[key])))) )

sys.exit()
