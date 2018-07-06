#!/usr/bin/env python
"""
collapsing.py

Function to perform collapsing analysis on simulated data. Gets disease status from .smp,
genotype data from .spl, returns FET p-value
"""

from scipy import stats

def collapse(smp_file, spl_file):
    # create sample dict samps[sample_name] = [disease status, qual_staus]
    samps = {}
    with open(smp_file, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            samps[line[0]] = [line[1], 0]
    # get qual status for each individual
    # later can be adapted to get genotypes if weighting scheme, etc
    # can also be adapted for recessive analysis
    with open(spl_file, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            indiv = line[0]
            if indiv in samps:
                if len(line) > 1:
                    samps[indiv][1] += 1
    qcase, qctrl, uqcase, uqctrl = 0, 0, 0, 0
    for k, v in samps.iteritems():
        if v[0] == '0':
            if v[1] == 0:
                uqctrl += 1
            else:
                qctrl += 1
        elif v[0] == '1':
            if v[1] == 0:
                uqcase += 1
            else:
                qcase += 1
    pval = stats.fisher_exact([[qcase, uqcase],[qctrl, uqctrl]])[1]
    #return (pval, smp_file, spl_file, qcase, uqcase, qctrl, uqctrl)
    return pval
