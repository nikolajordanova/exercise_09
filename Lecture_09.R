rm(list=ls()) #clean, clc, close all
# ***************
# R version 4.4.2 / win
# author: Nikola Jordánová
# *************

# Path
setwd('V:/MPA-PRG/exercise_09') # set working directory

library(Biostrings)


Overlap <- function(sequence1, sequence2){
  min_length <- min(nchar(sequence1), nchar(sequence2))
  max_overlap <- 0
  for (i in min_length:1){
    end_seq1 <- substr(sequence1, nchar(sequence1) - i + 1, nchar(sequence1))
    start_seq2 <- substr(sequence2, 1, i)
    if (end_seq1 == start_seq2){
      max_overlap <- i
      break
    }
  }
  return (max_overlap)
}

Overlap("AGACCTGCCG", "GCCGGAATAC")


OverlapMatrix <- function(S){
  len_S <- length(S)
  matrix_S <- matrix(0, len_S, len_S)
  
  for (i in 1:len_S){
    for (j in 1:len_S){
      if (i != j){
        matrix_S[i,j] <- Overlap(S[[i]], S[[j]])
      }
    }
  }
  return (matrix_S)
  
}


S <- DNAStringSet(c("CATGC", "CTAAGT", "GCTA", "TTCA", "ATGCATC"))
OverlapMatrix(S)


GreedySuperstring <-function(S){
  # S - DNAStringSet object of strings (reads)
  
  while (length(S) > 1){
    overlapMat <- OverlapMatrix(S)
    max_overlap <- max(overlapMat)
    if (max_overlap == 0){
      return (S)
    }
    else{
      position_of_max <- which(overlapMat == max_overlap, arr.ind = TRUE)[1, ]
      i <- position_of_max[1]
      j <- position_of_max[2]
      
      seq1 <- as.character(S[[i]])
      seq2 <- as.character(S[[j]])
      
      merged_sequences <- paste0(seq1, substr(seq2, max_overlap + 1, nchar(seq2))) # z druhé sekvence pouze část za překryvem
      
      S <- S[-c(i,j)]

      S <- append(S, DNAStringSet(merged_sequences))
      
      
    }
    
  }
  return (S)
}
    

# připravíme DNA sekvence
S <- DNAStringSet(c(
  "CATGC", "CTAAGT", "GCTA", "TTCA", "ATGCATC"
))

# zavolání funkce
GreedySuperstring(S)



