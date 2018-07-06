#!/usr/bin/env python
"""
collapsing.py

Function to perform collapsing analysis on simulated data. Gets disease status from .smp,
genotype data from .spl, returns FET p-value

This is an optimized version of the original
"""

from scipy import stats
import numpy as np
import random

def calc_p(burdens, phenos, lookup, pairs):
    case_indices = np.where(phenos == 1)[0] 
    ctrl_indices = np.where(phenos == 0)[0]
    qcase = np.where(burdens[case_indices] > 0)[0].shape[0]
    qctrl = np.where(burdens[ctrl_indices] > 0)[0].shape[0]
    uqcase = len(case_indices) - qcase
    uqctrl = len(ctrl_indices) - qctrl
    if (qcase, qctrl) in pairs:                                        
        pval = lookup[qcase, qctrl]                                    
    else:                                                              
        pval = stats.fisher_exact([[qcase, uqcase],[qctrl, uqctrl]])[1]
        lookup[qcase, qctrl] = pval                                    
        pairs.add((qcase, qctrl))                                      
    return float(pval), pairs                                          

def collapse(smp_file, spl_file, reg_file, nperms):
    loci = []
    # make dict of regions
    with open(reg_file, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            start, end = int(line[0]), int(line[1])
            loci.extend(range(start, end + 1))
    # create sample set (to check for inclusion) and sample list (for indexing)
    samp_set = set()
    samp_list = []
    # create phenotype vector
    phenos = []
    with open(smp_file, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            samp_set.add(line[0])
            samp_list.append(line[0])
            phenos.append(int(line[1]))
    phenos = np.array(phenos)
    ncase = np.where(phenos == 1)[0].shape[0]
    nctrl = np.where(phenos == 0)[0].shape[0]
    gt_matrix = np.zeros((len(samp_set), len(loci)))
    qual_case = 0
    qual_ctrl = 0
    # get gts for each individual
    with open(spl_file, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            indiv = line[0]
            if indiv in samp_set:
                if len(line) == 1:
                    pass
                else:
                    sample_index = samp_list.index(indiv)
                    phenotype = phenos[sample_index]
                    if phenotype == 1:
                        qual_case += 1
                    else:             
                        qual_ctrl += 1
                    for locus in line[1:]:
                        locus = int(locus)
                        locus_index = loci.index(locus)
                        gt_matrix[sample_index,locus_index] += 1
    #results = [None] * (nperms + 1)
    qual_sum = qual_case + qual_ctrl
    max_case = min([qual_sum, ncase])
    max_ctrl = min([qual_sum, nctrl])
    lookup = np.ones((max_case + 1, max_ctrl + 1))
    pairs = set()
    burdens = np.sum(gt_matrix, axis=1)
    prefix = spl_file[:-3]
    p_file = '{prefix}collapsing_pvals.txt'.format(prefix=prefix)
    with open(p_file, 'w') as outfile:
        obs_p, pairs = calc_p(burdens, phenos, lookup, pairs)
        #results[0] = obs_p
        outfile.write(str(obs_p) + '\n')
        for i in range(nperms):
            perm_phenos = random.sample(phenos, len(phenos))
            exp_p, pairs = calc_p(burdens, np.array(perm_phenos), lookup, pairs)
            #results[i+1] = exp_p
            outfile.write(str(exp_p) + '\n')
    #return results
    return "done"
