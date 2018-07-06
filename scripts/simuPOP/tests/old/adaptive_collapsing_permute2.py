#!/usr/bin/env python
"""
collapsing.py

Function to perform collapsing analysis on simulated data. Gets disease status from .smp,
genotype data from .spl, returns FET p-value

This is an optimized version of the original (significantly faster)
"""

from scipy import stats
import numpy as np
import random

def calc_p(burdens, phenos):
    case_indices = np.where(phenos == 1.0)[0] 
    ctrl_indices = np.where(phenos == 0.0)[0]
    qcase = np.where(burdens[case_indices] > 0)[0].shape[0]
    qctrl = np.where(burdens[ctrl_indices] > 0)[0].shape[0]
    uqcase = len(case_indices) - qcase
    uqctrl = len(ctrl_indices) - qctrl
    pval = stats.fisher_exact([[qcase, uqcase],[qctrl, uqctrl]])[1]
    return float(pval)

def collapse(smp_file, spl_file, reg_file, nperms):
    # make dict of regions
    regions = {}
    with open(reg_file, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            coeff = float('.' + line[2].strip('sel'))
            start, end = int(line[0]), int(line[1])
            regions[coeff] = [start, end]
    # sorted selection coefficients from most to least intolerant
    coeffs = sorted(regions.iterkeys(), reverse=True)
    # determine gene length for dimensions of GT matrix and make list of loci for indexing
    loci = []
    for k, v in regions.iteritems():
        reg_loci = range(v[0], v[1] + 1)
        loci.extend(reg_loci)
    loci = sorted(loci)
    # reset the positions based on relative position
    for k, v in regions.iteritems():
        v[0] = loci.index(v[0])
        v[1] = loci.index(v[1])
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
    results = [[None] * len(coeffs) for x in range(nperms + 1)]
    # the adaptive part - start with most intolerant region, then add on regions
    for i in range(len(coeffs)): 
        coeff_to_include = coeffs[0:i+1]
        sites_to_include = []
        # get list of sites in the regions that we want to include
        for coeff in coeff_to_include:
            sites = range(regions[coeff][0], regions[coeff][1] + 1)
            sites_to_include.extend(sites)
        sites_to_include = np.array(sites_to_include)
        gts = gt_matrix[:,sites_to_include]
        burdens = np.sum(gts, axis=1)
        obs_p = calc_p(burdens, np.array(phenos))
        results[0][i] = obs_p
        for j in range(nperms):
            perm_phenos = random.sample(phenos, len(phenos))
            exp_p = calc_p(burdens, np.array(perm_phenos))
            results[j+1][i] = exp_p
    min_results = [min(x) for x in results]
    return min_results
