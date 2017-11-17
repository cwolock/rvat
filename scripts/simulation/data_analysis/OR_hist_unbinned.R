#!/usr/bin/env Rscript

library(ggplot2)

#setwd('U:/assoc_test/simuPOP/schema2/uniform')
#ors <- read.table('regions.OR_dist.txt', header=TRUE, sep = '\t', as.is = TRUE)

# parse input parameters
args <- commandArgs(T)
in.file <- args[1]
out.path <- args[2]

ors <- read.table(in.file, header=TRUE, sep = '\t', as.is = TRUE)

pdf(out.path)
p <- ggplot(ors, aes(OR)) + 
  geom_histogram(binwidth=2) +
  ggtitle('Distribution of odds ratios for causal variants') +
  xlab('Odds ratio') + 
  ylab('Frequency') #+
  #ylim(0.00,1.00)
p
invisible(dev.off())
