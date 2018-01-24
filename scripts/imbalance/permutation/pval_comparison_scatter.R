#!/usr/bin/env Rscript

library(ggplot2)

setwd('U:/assoc_test/imbalance/permutation_based/comparison/balanced100')
matched <- read.table('matched_pvals.txt', header=TRUE, sep = '\t', as.is = TRUE)
matched <- matched[(matched$fet < 0.05 | matched$sum < 0.05),]

p <- ggplot(matched, aes(fet, sum)) + 
  geom_point() + 
  xlab('FET') + 
  ylab('Sum')
p

cor(matched$fet, matched$sum)
