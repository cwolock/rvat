#!/usr/bin/env python
"""
collapsing.py

Function to perform collapsing analysis on simulated data. Gets disease status from .smp,
genotype data from .spl, returns FET p-value

This version uses a lookup table
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
    #if q_case + 1 <= lookup.shape[0] and q_ctrl + 1 <= lookup.shape[1]:
    if (qcase, qctrl) in pairs:
        pval = lookup[qcase, qctrl] 
    else:
        pval = stats.fisher_exact([[qcase, uqcase],[qctrl, uqctrl]])[1]
        lookup[qcase, qctrl] = pval
        pairs.add((qcase, qctrl))
    return float(pval), pairs

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
    phenos = np.array(phenos)
    ncase = np.where(phenos == 1)[0].shape[0]
    nctrl = np.where(phenos == 0)[0].shape[0]
    gt_matrix = np.zeros((len(samp_set), len(loci)))
    # for lookup table
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
                    # calculate observed qual cases and ctrls for lookup table
                    phenotype = phenos[sample_index]
                    if phenotype == 1:
                        qual_case += 1
                    else: 
                        qual_ctrl += 1
                    for locus in line[1:]:
                        locus = int(locus)
                        locus_index = loci.index(locus)
                        gt_matrix[sample_index,locus_index] += 1
    #results = [[None] * len(coeffs) for x in range(nperms + 1)]
    # maximum number of qualified cases and controls (given permutation) in this gene
    qual_sum = qual_case + qual_ctrl
    max_case = min([qual_sum, ncase])
    max_ctrl = min([qual_sum, nctrl])
    lookup = np.ones((max_case + 1, max_ctrl + 1))
    # make list of (qual case, qual ctrl) tuples for FET lookup table
    #pairs = []
    pairs = set()
    #for i in xrange(max_case + 1):
    #    for j in xrange(max_ctrl + 1):
    #        pairs.append((i, j))
    #for pair in pairs:
    #    odds, p = stats.fisher_exact([[pair[0], ncase-pair[0]],[pair[1],nctrl-pair[1]]])
    #    lookup[pair[0],pair[1]] = p
    
    # the adaptive part - start with most intolerant region, then add on regions
    prefix = spl_file[:-3]
    p_file = '{prefix}adaptive_pvals.txt'.format(prefix=prefix)
    with open(p_file, 'w') as outfile:
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
            obs_p, pairs = calc_p(burdens, phenos, lookup, pairs)
            #results[0][i] = obs_p
            outfile.write(str(obs_p) + '\n')
            for j in range(nperms):
                perm_phenos = random.sample(phenos, len(phenos))
                exp_p, pairs = calc_p(burdens, np.array(perm_phenos), lookup, pairs)
                #results[j+1][i] = exp_p
                outfile.write(str(exp_p) + '\n')
        #min_results = [min(x) for x in results]
    #return min_results
    return "done"
