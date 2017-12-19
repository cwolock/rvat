#!/usr/bin/env Rscript

library(ggplot2)

# parse input parameters
args <- commandArgs(T)
in.file <- args[1]
out.path <- args[2]

vpk <- read.table(in.file, header=TRUE, sep = '\t', as.is = TRUE)

#OR$bin <- factor(as.character(OR$bin), 
#      levels=ordering)

pdf(out.path)
p <- ggplot(vpk, aes(vpk)) + 
  geom_histogram(binwidth=0.5) +
  ggtitle('Distribution of segregating sites per kb') +
  xlab('Seg. sites per kb') + 
  ylab('Frequency') #+
  #xlim(115, 200)
  #ylim(0.00,1.00)
p
invisible(dev.off())
