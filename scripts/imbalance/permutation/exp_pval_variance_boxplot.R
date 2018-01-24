#!/usr/bin/env Rscript

library(ggplot2)

setwd('U:/assoc_test/imbalance/permutation_based/comparison/balanced100')
coll_vars <- read.table('maf.001_delta.1.coll.variance.txt', header=FALSE, sep = '\t', as.is = TRUE)
colnames(coll_vars) <- 'variance'
coll_vars$test <- 'coll'
sum_vars <- read.table('maf.001_delta.1.sum.variance.txt', header=FALSE, sep='\t', as.is = TRUE)
colnames(sum_vars) <- 'variance'
sum_vars$test <- 'sum'
variances <- rbind(coll_vars, sum_vars)
#variances <- variances[variances$variance > 0.00002,]

p <- ggplot(variances, aes(test, variance)) + 
  geom_boxplot() +
  ggtitle('Distribution of odds ratios for causal variants') +
  xlab('Region') + 
  ylab('Odds ratios') #+
#ylim(0.00,1.00)
p
