#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# setwd('/home/henri/Repos/basf/dash_dashboard')
# test if there is at least one argument: if not, return a error
if (length(args)!=3) {
  stop("Pass in 1. filename, 2. sequences, 3. sequences.n", call.=FALSE)
}

library(DiffLogo)
library(Biostrings)

# pass in string of sequences, comma separated
sequences1 = args[2]
sequences2 = args[3]

# convert to character class
sequences1 <- strsplit(sequences1, ',')
sequences1 <- unlist(sequences1)

sequences2 <- strsplit(sequences2, ',')
sequences2 <- unlist(sequences2)

# make Position Weight Matrices
# from https://www.biostars.org/p/107356/
string_set1 <- AAStringSet(sequences1)
PSSM1 <- consensusMatrix(string_set1, as.prob = TRUE)
# colnames(PSSM1) = as.character(unlist(1:length(PSSM1[1,])))

string_set2 <- AAStringSet(sequences2)
PSSM2 <- consensusMatrix(string_set2, as.prob = TRUE)

# from the DiffLogo.R tutorial script
# actually want the ratio to adapt to window size...
widthToHeightRatio = nchar(sequences1[1])/8;
size = 2
resolution = 100
width = size * widthToHeightRatio
height = size * 2


# plot DiffSeq logos to file
png(file = args[1], res = resolution,
 width = width * resolution, height = height * resolution)

diffLogoFromPwm(PSSM1, PSSM2, alphabet=ASN)
graphics.off()

