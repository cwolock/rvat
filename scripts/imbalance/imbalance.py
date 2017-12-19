#!/usr/bin/env python
"""
imbalance_testing.py

Test for case/control enrichment under the null (no causal variants)
with imbalanced cohort size
"""

from scipy import stats
import argparse
import numpy as np
import random
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('darkgrid')

def run(gene_size, n_case, n_ctrl, delta, max_maf, iters, plot_file):
    pvals = []
    enriches = []
    for i in range(iters):
        gene = np.array([0.0] * gene_size)
        variants = random.sample(range(gene_size), int(np.ceil(delta*gene_size)))
        afs = stats.uniform.rvs(loc=0, scale=max_maf, size=len(variants))
        gene[variants] = afs
        cases = np.zeros((n_case, gene_size), dtype=np.float64)
        ctrls = np.zeros((n_ctrl, gene_size), dtype=np.float64)
        for pos, af in zip(variants, afs):
            case_gts = stats.bernoulli.rvs(p=af, size=n_case)
            ctrl_gts = stats.bernoulli.rvs(p=af, size=n_ctrl)
            cases[:,pos] = case_gts
            ctrls[:,pos] = ctrl_gts
        qcase = len(set(np.where(cases > 0)[0]))
        qctrl = len(set(np.where(ctrls > 0)[0]))
        print qcase, qctrl
        pval = stats.fisher_exact([[qcase, n_case-qcase],[qctrl, n_ctrl-qctrl]])[1]
        pvals.append(-np.log10(pval))
        enriches.append(qcase/float(qcase+n_case) - qctrl/float(qctrl+n_ctrl))
        
    fig = plt.figure(figsize=(12,12))
    plt.xlabel('-log10(p)', fontsize=20)
    plt.ylabel('Case enrichment', fontsize=20)
    dataAx = fig.add_subplot(1,1,1)
    dataAx.plot(pvals, enriches, 'r.', label='_nolegend_', markersize=10)
    plt.tight_layout()
    plt.savefig(plot_file)

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
    parser.add_argument('--plot_file', help='Specify name of plot file')
    args = parser.parse_args()
    run(args.gene_size, args.n_case, args.n_ctrl, args.delta, 
        args.max_maf, args.iters, args.plot_file)
