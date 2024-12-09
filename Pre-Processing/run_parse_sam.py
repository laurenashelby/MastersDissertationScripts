#!/share/apps/python/3.8.6/intel/bin/python

"""
Title:        SAM to tsv
Author:       Ziyue Cheng
"""

###########################################
# This script is used to parse SAM file
# 
# INPUT: SAM
# OUTPUT: tsv with 6 columns
#      1. STRAND: +/- (if NH:1, FLAG can only be 0 or 16)
#      2. RNAME: name of reference sequence/chrom
#      3. POS5: position of 5' end
#      4. CIGAR: CIGAR from SAM file
#      5. POS3: position of 3' end
#      6. LEN: length of the read
#
# NOTE: you may want to change the first line to your system python interpreter
###########################################

import sys

NUMBER_CHAR = '0123456789'

def get_ref_length(cigar:str):
    """
    Get the match length on reference based on CIGAR column of SAM file
    """
    length = 0
    number = ""
    for l in cigar:
        if l in NUMBER_CHAR:
            number += l
        else:
            # Common CIGAR character: M I D N
            # I means insertion to reference
            # M match, D deletion, N skip (intron in RNA seq)
            # so to get the match length on reference, we should ignore insertion number
            if l != "I":
                length += int(number)
            number = ""
    return length

split_char = None

sys.stdout.write("STRAND\tRNAME\tPOS5\tCIGAR\tPOS3\tLEN\n")

for line in sys.stdin:
    # each line in std input
    if not split_char: # if split_char is None
        # if table in line, then the separater is table symbol
        split_char = "\t" if "\t" in line else " "

    elements = line.split(split_char)
    
    sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
        "+" if elements[1] == "0" else "-", # 0 plus strand, 16 negative strand 
        elements[2], # ref name
        elements[3], # 5' position
        elements[5], # CIGAR
        int(elements[3]) + get_ref_length(elements[5]) - 1, # 3' position
        len(elements[9]) # length of read
    ))

