#!/usr/bin/env python
"""
weighted_sum.py

Perform Sum test weighted by MAF
"""

from scipy import stats
from utilities import *
import argparse
import numpy as np

def run(gene_size, n_case, n_ctrl, delta, max_maf, iters, output_prefix):
    qq_file = output_prefix + '.qq_plot.pdf'
    log_file = output_prefix + '.log'
    write_log(log_file, gene_size, n_case, n_ctrl, delta, max_maf)
    pvals = []
    for i in range(iters):
        cases, ctrls, variants, afs = generate_null_gts(gene_size, n_case, 
                                                        n_ctrl, delta, max_maf)
        for pos, af in zip(variants, afs):
            cases[:,pos] = cases[:,pos]/float(af)
            ctrls[:,pos] = ctrls[:,pos]/float(af)
        case_sums = np.sum(cases, 1)
        ctrl_sums = np.sum(ctrls, 1)
        X_vector = np.append(case_sums, ctrl_sums)
        Y_vector = np.array([1.0]*n_case + [0.0]*n_ctrl)
        X_bar = np.mean(X_vector)
        Y_bar = n_case/float(n_case + n_ctrl)
        Y_diff = Y_vector - Y_bar
        X_diff = X_vector - X_bar
        U = sum(np.multiply(Y_diff, X_vector))
        V = Y_bar * (1-Y_bar) * sum(np.multiply(X_diff, X_diff))
        score = U * 1/float(V) * U
        pval = 1 - stats.chi2.cdf(score, 1)
        if str(pval) == 'nan': pval = 1
        pvals.append(pval)
    #pvals = [x for x in pvals if str(x) != 'nan']
    pvals = sorted(pvals)
    #for p in pvals:
    #    if p == 'nan': print 
    plot_qq(pvals, qq_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--gene_size', type=int, help='Specify gene length')
    parser.add_argument('--n_case', type=int, help='Specify number of cases')
    parser.add_argument('--n_ctrl', type=int, help='Specify number of controls')
    parser.add_argument('--delta', type=float, help='Specify % of sites with a variant')
    parser.add_argument('--max_maf', type=float, help='Specify max MAF')
    parser.add_argument('--iters', type=int, help='Specify number of iterations')
    parser.add_argument('--output_prefix', help='Specify output file prefix')
    args = parser.parse_args()
    run(args.gene_size, args.n_case, args.n_ctrl, args.delta, 
        args.max_maf, args.iters, args.output_prefix)
