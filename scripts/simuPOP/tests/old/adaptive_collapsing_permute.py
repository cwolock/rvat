#!/usr/bin/env python
"""
collapsing.py

Function to perform collapsing analysis on simulated data. Gets disease status from .smp,
genotype data from .spl, returns FET p-value
"""

from scipy import stats
import numpy as np
import random

def calc_p(coeffs, regions, loci, phenos, gt_matrix):
    results = []
    # the adaptive part - start with most intolerant region, then add on regions
    for i in range(len(coeffs)): 
        coeff_to_include = coeffs[0:i+1]
        sites_to_include = []
        # get list of sites in the regions that we want to include
        for coeff in coeff_to_include:
            sites = range(regions[coeff][0], regions[coeff][1] + 1)
            sites_to_include.extend(sites)
        sites_to_include = set(sites_to_include)
        site_indices = []
        # get site (column) indices of sites to include
        for site in sites_to_include:
            site_index = loci.index(site)
            site_indices.append(site_index)
        site_indices = np.array(site_indices)
        # track qualified cases and controls
        qcase, qctrl, uqcase, uqctrl = 0, 0, 0, 0
        for i in range(gt_matrix.shape[0]):
            gts = gt_matrix[i,site_indices]
            if sum(gts) > 0:
                if phenos[i] == 0:
                    qcase += 1
                elif phenos[i] == 1:
                    qctrl += 1
            else:
                if phenos[i] == 0:
                    uqctrl += 1
                elif phenos[i] == 1:
                    uqcase += 1
        pval = stats.fisher_exact([[qcase, uqcase],[qctrl, uqctrl]])[1]
        results.append(pval)
    return min(results) 

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
    pval_file = spl_file[:-3] + 'pvals.txt'
    with open(pval_file, 'w') as outfile:
        obs_p = calc_p(coeffs, regions, loci, phenos, gt_matrix)
        outfile.write(str(obs_p) + '\n')
        #expected_dist = []
        for i in range(nperms):
            perm_phenos = random.sample(phenos, len(phenos))
            expected = calc_p(coeffs, regions, loci, perm_phenos, gt_matrix)
            outfile.write(str(expected) + '\n')
            #expected_dist.append(expected)
        #all_p = [obs_p] + expected_dist
    #expected_dist = np.array(expected_dist)
    #perm_p = np.where(expected_dist < obs_p )[0].shape[0] / float(nperms)
    #return perm_p
    return 'success'
