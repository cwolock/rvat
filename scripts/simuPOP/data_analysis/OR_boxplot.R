#!/usr/bin/env Rscript

library(ggplot2)

#setwd('U:/assoc_test/simuPOP/schema2/uniform')
#ors <- read.table('regions.OR_dist.txt', header=TRUE, sep = '\t', as.is = TRUE)

# parse input parameters
args <- commandArgs(T)
in.file <- args[1]
out.path <- args[2]

ors <- read.table(in.file, header=TRUE, sep = '\t', as.is = TRUE)

ordering <- c('all_regions', 'reg1', 'reg2', 'reg3', 'reg4')
ors$regs <- factor(as.character(ors$reg), levels=ordering)

pdf(out.path)
p <- ggplot(ors, aes(reg, OR)) + 
  geom_boxplot() +
  ggtitle('Distribution of odds ratios for causal variants') +
  xlab('Region') + 
  ylab('Odds ratios') #+
  #ylim(0.00,1.00)
p
invisible(dev.off())
