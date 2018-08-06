# 08.03.2018
# Creator: Susana Wilson Hawke
# Identify variants in IDRs in 3' end of protein (last IDR in sequence)

library(dplyr)
library(data.table)
library(protr)
library(factoextra)
library(cluster) 
library(gridExtra)
library(plyr)
library(readr)
library(gridExtra)
library(grid)
library(WriteXLS)

#------------------------------------------------------------------------------#

#                             Define functions                                 #

#------------------------------------------------------------------------------#
getSequences <- function(x){
  
  sequences <- c()
  ensembl <- c()
  genenames <- c()
  
  for (i in 1:nrow(x)){
    if (i%%2 == 0){
      sequence = as.character(x$V1[i])
      sequences <- c(sequences, sequence)
    } else if(i%%2 == 1 ){
      id = as.character(x$V1[i])
      id <- strsplit(id, ">")[[1]][2]
      ensembl <- c(ensembl, id)
      name = as.character(x$V2[i])
      genenames <- c(genenames, name)
    }
  }
  
  return(data.frame(ensembl, genenames, sequences))
}

#------------------------------------------------------------------------------#

#                           Read and preprocess datasets                       #

#------------------------------------------------------------------------------#

# read command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least two argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  stop("At least two arguments must be supplied (input file).n", call.=FALSE)
} else {
  args[3] = "out.txt"
}

df <- read.csv(args[1], sep = "\t", col.names = c("peptide.id", "aa.pos1", "aa.pos2", 
                                                                           "genename", "locus", "disease", "peptide.id1",
                                                                           "idr.start", "idr.end", "misc", "genename2"))

setwd("~/Desktop/Young_Lab_2018/ClinVar_D2P2/starting_files/")
sequences <- read.csv("Ensembl_proteome.clean.fa", sep = "|", header = FALSE)
clean.table <- getSequences(sequences)
write.table(clean.table, file = "entire.proteome.fastafile", sep = "\t", quote = FALSE, row.names = FALSE)


#------------------------------------------------------------------------------#

#                           Read and clean datasets                            #

#------------------------------------------------------------------------------#

sequences <- read.csv("entire.proteome.fastafile", sep = "\t", header = TRUE)
sequences <- dplyr::rename(sequences, peptide.id = ensembl, genename = genenames)

aa.pos1 = rep("a", nrow(sequences))
aa.pos2 = rep("a", nrow(sequences))
locus = rep("a", nrow(sequences))
disease = rep("a", nrow(sequences))
peptide.id1 = rep("a", nrow(sequences))
idr.start = rep("a", nrow(sequences))
idr.end = rep("a", nrow(sequences))
misc = rep("a", nrow(sequences))
genename2 = rep("a", nrow(sequences))

sequences$aa.pos1 = aa.pos1
sequences$aa.pos2 = aa.pos2
sequences$locus = locus
sequences$disease = disease
sequences$peptide.id1 = peptide.id1
sequences$idr.start = idr.start
sequences$idr.end = idr.end
sequences$misc = misc
sequences$genename1 = genename2

temp <- left_join(df, sequences, by = "peptide.id")

peptide.id <- temp$peptide.id
aa.pos1 <- temp$aa.pos1.x
aa.pos2 <- temp$aa.pos2.x
genename <- temp$genename.x
locus <- temp$locus.x
disease <- temp$disease.x
peptide.id1 <- temp$peptide.id1.x
idr.start <- temp$idr.start.x
idr.end <- temp$idr.end.x
misc <- temp$misc.x
sequences <- temp$sequences

final <- data.frame(peptide.id, aa.pos1, aa.pos2, genename, locus, disease, peptide.id1, idr.start, idr.end, misc, sequences) 
df.seq <- final

df.seq <- mutate(df.seq, strsplit = strsplit(as.character(sequences), ""))
length <- unlist(lapply(df.seq$strsplit, function(x){length(unlist(x))}))

df.seq$length <- length
df.seq$halfway <- df.seq$length/2

df.seq <- df.seq %>%  subset(as.numeric(as.character(df.seq$aa.pos1)) > df.seq$halfway)
final <- subset(df.seq, select = -c(sequences,strsplit))

# create filename for final file
filename <- c(args[2],".3prime.xls" )
filename <- paste0(filename, sep = "", collapse = "")

# write an excel file with results
WriteXLS(final, ExcelFileName = filename, row.names = FALSE)




