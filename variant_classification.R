# 07.30.2018
# Creator: Susana Wilson Hawken

#-------------------------------------------------------------------#

#      analysis of STOPGAIN mutations outside of DNA binding        #
#      domain in transcription factors                              #

#-------------------------------------------------------------------#

#download dplyr package
library(dplyr)
library(data.table)
library(protr)
library(factoextra)
library(cluster) 
library(gridExtra)
library(plyr)
library(WriteXLS)

setwd("~/Desktop/Young_Lab_2018/ClinVar_D2P2/variant_files/inIDR/in3prime/NotInDomains")
d <- read.csv("../NONSYNONYMOUS.3prime.csv", sep = "\t", header = FALSE)
g <- read.csv("aa.genetable_clean.bed", sep = "\t", header = FALSE)

d <- as.data.table(d)
g <- as.data.table(g)

d_clean <- d[!duplicated(d),]
g_clean <- g[!duplicated(g),]

d_clean <- dplyr::rename(d_clean, protein.id = V1, aa.position = V2, aa.position2 = V3, genename = V4, locus = V5, disease = V6, protein.id2 = V7, idr.start= V8, idr.end = V9, misc = V10, length = V11, halfway = V12)
g_clean <- dplyr::rename(g_clean, protein.id = V1, aa.position = V2, aa.position2 = V3, genename = V4, sequence = V5)

locus <- rep("0", nrow(g_clean))
disease <- rep("a", nrow(g_clean))
protein.id2 <- rep("0", nrow(g_clean))
idr.start <- rep("0", nrow(g_clean))
idr.end <- rep("0", nrow(g_clean))
misc <- rep("0", nrow(g_clean))
length <- rep("0", nrow(g_clean))
halfway <- rep("0", nrow(g_clean))

g_clean$locus <- locus
g_clean$disease <- disease
g_clean$protein.id2 <- protein.id2
g_clean$idr.start <- idr.start 
g_clean$idr.end <- idr.end
g_clean$misc <- misc
g_clean$length <- length
g_clean$halfway <- halfway

sequence <- rep("a", nrow(d_clean))

d_clean$sequence = sequence

temp <- left_join(g_clean, d_clean, by = "protein.id")
temp <- temp[!is.na(temp$aa.position.y),]

protein.id <- temp$protein.id
aa.pos <- temp$aa.position.x
aa.pos2 <- temp$aa.position2.x
idr.start <- temp$idr.start.y
idr.end <- temp$idr.end.y
mut.pos <- temp$aa.position.y
genename <- temp$genename.x
length <- temp$length.y

f <- data.table(protein.id, aa.pos, aa.pos2, idr.start, idr.end, mut.pos, genename, length)
final <- f[which(as.numeric(as.character(aa.pos2)) < as.numeric(as.character(mut.pos))),]
write.table(final, file = "NONSYNONYMOUS.after.domain.txt", sep = "\t", quote = FALSE, row.names = FALSE)
WriteXLS(final, ExcelFileName = "NONSYNONYMOUS.after.domain.xls", row.names = FALSE)
