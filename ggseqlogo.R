#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# setwd('/home/henri/Repos/basf/dash_dashboard')
# test if there is at least one argument: if not, return a error
if (length(args)!=2) {
  stop("Pass in 1. filename, 2. sequence.n", call.=FALSE)
}

library(Biostrings)
require(ggplot2)
require(ggseqlogo)

# pass in string of sequences, comma separated
sequences1 = args[2]

# convert to character class
sequences1 <- strsplit(sequences1, ',')
sequences1 <- unlist(sequences1)

# from the DiffLogo.R tutorial script
# actually want the ratio to adapt to window size...
widthToHeightRatio = nchar(sequences1[1])/8;
size = 2
resolution = 100
width = size * widthToHeightRatio
height = size


# plot DiffSeq logos to file
png(file = args[1], res = resolution,
    width = width * resolution, height = height * resolution)

ggseqlogo( sequences1, seq_type='aa', method = 'bits' ) #, col_scheme = 'hydrophobicity' )

graphics.off()

