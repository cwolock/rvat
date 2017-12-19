#!/usr/bin/env Rscript

library(ggplot2)

# parse input parameters
args <- commandArgs(T)
OR.file <- args[1]
out.path <- args[2]

OR <- read.table(OR.file, header = TRUE, sep = '\t', as.is = TRUE)

ordering <- c('(1, 1.1)','(1.1, 1.2)','(1.2, 1.3)','(1.3, 1.4)','(1.4, 1.5)','(1.5, 1.6)','(1.6, 2)','(2, 3)','(3, 5)','(5, 10)','(10, 100)')

OR$bin <- factor(as.character(OR$bin), 
      levels=ordering)

pdf(out.path)
p <- ggplot(OR, aes(x=bin, y=frac, fill=type)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle('Odds ratio distributions for rare and common variants') +
  xlab('OR bin') + 
  ylab('Fraction of variants') +
  ylim(0.00,1.00)
p
invisible(dev.off())
