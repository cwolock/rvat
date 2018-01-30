#!/usr/bin/env python
"""
collapsing.py

Function to perform collapsing analysis on simulated data. Gets disease status from .smp,
genotype data from .spl, returns FET p-value
"""

from scipy import stats

def collapse(smp_file, spl_file, reg_file):
    # make dict of regions
    regions = {}
    with open(reg_file, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            coeff = float('.' + line[2].strip('sel'))
            print coeff
            start, end = int(line[0]), int(line[1])
            regions[coeff] = [start, end]
    # sorted selection coefficients from most to least intolerant
    coeffs = sorted(regions.iterkeys(), reverse=True)
    # create sample dict samps[sample_name] = [disease status, gts]
    samps = {}
    with open(smp_file, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            samps[line[0]] = [line[1], []]
    # get gts for each individual
    with open(spl_file, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            indiv = line[0]
            if indiv in samps:
                if len(line) == 1:
                    pass
                else:
                    samps[indiv][1].append(map(int, line[1:]))
    results = []
    # the adaptive part - start with most intolerant region, then add on regions
    for i in range(len(coeffs)): 
        coeff_to_include = coeffs[0:i+1]
        sites_to_include = []
        # get list of sites in the regions that we want to include
        for coeff in coeff_to_include:
            sites = range(regions[coeff][0], regions[coeff][1] + 1)
            sites_to_include.append(sites)
        # track qualified cases and controls
        qcase, qctrl, uqcase, uqctrl = 0, 0, 0, 0
        for k, v in samps.iteritems():
            qvs = [x for x in v[1] if x in sites_to_include]
            if v[0] == '0':
                if len(qvs) > 0:
                    qctrl += 1
                else:
                    uqctrl += 1  
            elif v[0] == '1':
                if len(qvs) > 0:
                    qcase += 1
                else:
                    uqcase += 1
        pval = stats.fisher_exact([[qcase, uqcase],[qctrl, uqctrl]])[1]
        results.append((coeff_to_include, pval))
    return results
