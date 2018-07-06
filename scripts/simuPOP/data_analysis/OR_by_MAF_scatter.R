#!/usr/bin/env Rscript

library(ggplot2)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#setwd('U:/assoc_test/simuPOP/schema2/linear/paf0.01')
#ors <- read.table('multi_region.OR_by_MAF.txt', header=TRUE, sep = '\t', as.is = TRUE)

# parse input parameters
args <- commandArgs(T)
in.file <- args[1]
out.path <- args[2]

ors <- read.table(in.file, header=TRUE, sep = '\t', as.is = TRUE)
ors$MAF <- log10(ors$MAF)
#ordering <- c('reg1', 'reg2', 'reg3', 'reg4')
#ordering <- c('reg4', 'reg3', 'reg2', 'reg1')
ordering <- c('sel01', 'sel001', 'sel0')
ors$regs <- factor(as.character(ors$reg), levels=ordering)

pdf(out.path)
p <- ggplot(ors, aes(MAF, OR)) + 
  geom_jitter(aes(color = reg, alpha=reg),width=0.1, height=2,shape=19) + #alpha=0.2
  #geom_jitter(data=subset(ors,reg=='0'), aes(x=MAF,y=OR,color=reg, alpha=reg), width=0.1,height=2,shape=19) +
  #scale_color_manual(values=c('grey85', 'grey60', 'grey35', '')) +
  #scale_color_manual(values=cbbPalette,name='Sel. Coeff.', labels = c('0', '0.001', '0.01', '0.1')) + 
  scale_color_manual(values=cbbPalette,name='Sel. Coeff.', labels = c('0', '0.001', '0.01')) + 
  #scale_alpha_manual(guide='none',values=c('0'=0.2,'0.001'=0.2,'0.01'=0.2,'0.1'=1)) + 
  #scale_shape(solid=FALSE) +
  ggtitle('Distribution of odds ratios for causal variants') +
  #scale_color_discrete(values=cbbPalette, name='Selection Coeff.', labels = c('0', '0.001', '0.01', '0.1')) + 
  xlab('log10(MAF)') + 
  ylab('Odds ratio') +
  #scale_alpha_manual(guide='none',values=c(1,0.5, 0.5, 0.5))#+
  scale_alpha_manual(guide='none',values=c(0.5,0.5, 1))#+
  #ylim(0.00,1.00)
p
invisible(dev.off())
