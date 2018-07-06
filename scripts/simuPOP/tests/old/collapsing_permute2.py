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

def calc_p(burdens, phenos):
    case_indices = np.where(phenos == 1)[0] 
    ctrl_indices = np.where(phenos == 0)[0]
    qcase = np.where(burdens[case_indices] > 0)[0].shape[0]
    qctrl = np.where(burdens[ctrl_indices] > 0)[0].shape[0]
    uqcase = len(case_indices) - qcase
    uqctrl = len(ctrl_indices) - qctrl
    pval = stats.fisher_exact([[qcase, uqcase],[qctrl, uqctrl]])[1]
    return pval

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
    gt_matrix = np.zeros((len(samp_set), len(loci)))
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
                    for locus in line[1:]:
                        locus = int(locus)
                        locus_index = loci.index(locus)
                        gt_matrix[sample_index,locus_index] += 1
    results = [None] * (nperms + 1)
    burdens = np.sum(gt_matrix, axis=1)
    obs_p = calc_p(burdens, np.array(phenos))
    results[0] = obs_p
    for i in range(nperms):
        perm_phenos = random.sample(phenos, len(phenos))
        exp_p = calc_p(burdens, np.array(perm_phenos))
        results[i+1] = exp_p
    return results
