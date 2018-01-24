#!/usr/bin/env python
"""
standard_collapsing.py

Perform standard collapsing under the null (no causal variants)
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
        qcase = len(set(np.where(cases > 0)[0]))
        qctrl = len(set(np.where(ctrls > 0)[0]))
        OR, pval = stats.fisher_exact([[qcase, n_case-qcase],[qctrl, n_ctrl-qctrl]])
        pvals.append(pval)
    pvals = sorted(pvals)
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
