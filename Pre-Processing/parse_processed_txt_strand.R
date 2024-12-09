
# =============================================================================
# Title:        Alignment of txt to Annotation with transcript info
# Author:       Ziyue Cheng
# =============================================================================

# configuration ----

#########################################
#
# INPUT: star/*_unique.txt.gz
# OUTPUT: annotation with read information
#
# This is a memory-intensive script.
# Please allocate at least more than 5GB of memory per core. Ideally, 8GB per core.
# And it takes a lot of time to run. Set it aside and have a cup of coffee.
#
#########################################

sample_names = c( #TOCHANGE
  "SRR10416856",  # 80S WT
  "SRR11553468",  # 40S WT
)

# optimize for SLURM on HPC
nCPU = Sys.getenv("SLURM_CPUS_PER_TASK")
# if $SLURM_CPUS_PER_TASK == None, thread_number = 4. Else thread_number = $SLURM_CPUS_PER_TASK
thread_number = if(nCPU == "") 4 else as.numeric(nCPU)

print(paste0("[", date(), "] START with ", thread_number, " cores"))
# Preparing data ----
library(tidyverse)
library(multidplyr)

# This function will return the read info
# 
# FORMAT
#       Read Length:Location    separated by comma
# e.g.  31x375,31x480,35x512,54x275
get_read_info = function(current_strand, current_txStart, current_txEnd, current_exonStarts, current_exonEnds, current_chrom) {
  # select reads in this transcript
  sample_current = sample %>%
    # NOTE: if it is cDNA, can we ignore the strand?
    filter(POS5 - 1 >= current_txStart, POS3 <= current_txEnd, RNAME == current_chrom)
  # filter(RNAME == current_chrom, STRAND == current_strand, POS5 - 1 >= current_txStart, POS3 <= current_txEnd)
  
  # If no reads, return NA
  if (nrow(sample_current) == 0)
    return(NA)
  
  # add transcript strand to sample_current table
  sample_current["txStrand"] = current_strand
  
  # count reads distance to 5'
  sample_current = sample_current %>%
    mutate(distance = pmap_int(
      # list(STRAND, POS5, POS3, LEN),
      list(txStrand, POS5, POS3, LEN), # use transcript strand not the read strand
      function(strand, readStart, readEnd, readLen) {
        # txStart + 1 is the transcription start in GENCODE, the readStart should minus 1 here.
        readStart = readStart - 1
        
        exonStarts = unlist(current_exonStarts)
        exonEnds   = unlist(current_exonEnds)
        exonCount  = length(exonStarts)
        
        readStartRank = ceiling(rank(c(readStart, exonStarts))[1])
        readEndRank   = floor(rank(c(readEnd, exonEnds))[1])
        exonLengths  = exonEnds - exonStarts
        
        # ENSURE the read is not in the intron
        # invalid 5' position
        # if (readStartRank == 0 || readStart > exonEnds[readStartRank-1]) {
        if (readStart > exonEnds[readStartRank-1]) {
          return(NA)
        }
        # invalid 3' position
        # if (readEndRank == exonCount + 1 || readEnd < exonStarts[readEndRank]) {
        if (readEnd < exonStarts[readEndRank]) {
          return(NA)
        }
        
        if (strand == "+") {
          # ---- ---- ---== ====
          #           |||          ||| is where we calculate
          distance_to_5 = readStart - exonStarts[readStartRank-1] + 1 # remember to plus 1
          
          if (readStartRank > 2) {
            # ---- ---- ---== ====
            # |||| ||||            ||| is where we calculate here
            distance_to_5 = distance_to_5 + sum(exonLengths[1:(readStartRank-2)]) 
          }
        } else {
          # ==== ==--- ---- ----
          #        |||           ||| is where we calculate
          distance_to_5 = exonEnds[readEndRank] - readEnd + 1 # remember to plus 1
          
          if(readEndRank < exonCount) {
            # ==== ==--- ---- ----
            #            |||| ||||  ||| is where we calculate here
            distance_to_5 = distance_to_5 + sum(exonLengths[(readEndRank+1):exonCount])
          }
        }

        distance_to_5
      }
    ))
  
  # filter out NA
  sample_current = sample_current %>%
    filter(!is.na(distance)) %>%
    select(LEN,distance) %>%
    arrange(LEN, distance)
  
  # If no reads, return NA
  if (nrow(sample_current) == 0)
    return(NA)
  
  # save the read length and distance to the
  # FORMAT: readLen x distance to 5' , seperated by comma
  # e.g.    31x375,31x480,35x512,54x275
  paste(sample_current$LEN,sample_current$distance, sep = "x", collapse = ",")
}

# do in a for loop ----
for(s in sample_names){
  print(paste0("[", date(), "] ", "Start analyzing ", s))
  # Read annotation ----
  stst_annotation = read.delim("APPRIS_ENSEMBL109_stst_uorf_u9-15_gene_transcript.tsv")
  
  # convert comma separated columns
  stst_annotation = stst_annotation %>%
    # split comma ,
    mutate(across(matches("s$"), ~ strsplit(., ","))) %>%
    # convert to int vector
    mutate(across(matches("(start|end|2CDS)s$"), ~ map(., as.numeric) ))
  
  # read processed file ----
  sample_path = paste0("star/", s,"_unique.txt.gz") #TOCHANGE
  output_path = paste0("parse-sam-double/", s,"_with_read_info_final.tsv.gz") #TOCHANGE
  sample = read.delim(sample_path, header = T, sep="\t")
  print(object.size(sample), units = "MB")
  
  # preparing for multithread ----
  cluster = new_cluster(thread_number)
  cluster_library(cluster, "tidyverse") # load library on clusters
  cluster_copy(cluster, "sample")
  cluster_copy(cluster, "get_read_info") # copy data and functions to clusters
  
  # add read info to this table ----
  stst_annotation = stst_annotation %>%
    partition(cluster) %>% # split data to clusters
    mutate(readInfo = pmap_chr(
      list(strand, txStart, txEnd, exonStarts, exonEnds, chrom),
      get_read_info
    )) %>%
    collect() # collect from clusters
  
  # save to file ----
  print(paste0("[", date(), "] ", "Start saving to ", output_path))
  stst_annotation = stst_annotation %>%
    mutate(across(matches("s$"), ~ map_chr(.,  ~ paste(., collapse=","))))
  
  write.table(stst_annotation, gzfile(output_path), row.names = F, quote = F, sep = "\t")

  # remove the cluster
  rm(cluster)
  gc()
}

print(paste0("[", date(), "] COMPLETE"))

