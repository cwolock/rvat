#/usr/bin/env python
"""
utilities.py

Helpful functions for imbalance analyses
"""

from scipy import stats
import numpy as np
import random
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('darkgrid')

def generate_null_gts(gene_size, n_case, n_ctrl, delta, max_maf):
    gene = np.array([0.0] * gene_size)
    #variants = random.sample(range(gene_size), int(np.ceil(delta*gene_size)))
    n_vars = stats.binom.rvs(n=gene_size, p=delta, loc=0, size=1)[0]
    variants = random.sample(range(gene_size), n_vars)
    afs = stats.uniform.rvs(loc=0, scale=max_maf, size=len(variants))
    gene[variants] = afs
    cases = np.zeros((n_case, gene_size), dtype=np.float64)
    ctrls = np.zeros((n_ctrl, gene_size), dtype=np.float64)
    for pos, af in zip(variants, afs):
        case_gts = stats.bernoulli.rvs(p=af, size=n_case)
        ctrl_gts = stats.bernoulli.rvs(p=af, size=n_ctrl)
        cases[:,pos] = case_gts
        ctrls[:,pos] = ctrl_gts
    return (cases, ctrls, variants, afs)

def plot_qq(pvals, qq_file):
    # generate uniform expected
    expected = sorted(stats.uniform.rvs(loc=0, scale=1, size=len(pvals)))
    #for a, b in zip(pvals, expected):
    #    print a,b
    #expected = np.arange(0,1,1/float(len(pvals)))
    
    # transform variables
    obs = -np.log10(pvals)
    exp = -np.log10(expected)
    
    # set limits of qq plot axes                                                   
    axisMax_obs = np.ceil(max(obs))                                                
    axisMax_exp = np.ceil(max(exp))                                                
                                                                                   
    # initialize qqplot with axes and labels                                       
    fig = plt.figure(figsize=(12,12))                                              
    plt.xlim([0, axisMax_exp])                                                     
    plt.xlabel('Expected -log10(p)', fontsize=20)                                  
    plt.ylim([0, axisMax_obs])                                                     
    plt.ylabel('Observed -log10(p)', fontsize=20)                                  
    plt.title('QQ Plot: Observed vs. expected p-values', fontsize=20)
    
    # change size of axis tick labels                                              
    plt.tick_params(axis='both', which='major', labelsize=12)                      
                                                                                   
    # plot the points, exp on x axis, obs on y axis                                
    dataAx = fig.add_subplot(1,1,1)                                                
    dataAx.plot(exp, obs, 'r.', label='_nolegend_', markersize=12)                 
                                                                                   
    # plot a diagonal line for comparison                                          
    lineAx = fig.add_subplot(1,1,1)                                                
    lineAx.plot([0,max(axisMax_obs,axisMax_exp)], [0,max(axisMax_obs,axisMax_exp)],
                 'b-', label='_nolegend_')

    plt.tight_layout()
    plt.savefig(qq_file)

def write_log(log_file, gene_size, n_case, n_ctrl, delta, max_maf):
    with open(log_file, 'w') as outfile:                                                  
        outfile.write('gene_size={}\nn_case={}\nn_ctrl={}\ndelta={}\nmax_maf={}\n'.format(
            gene_size, n_case, n_ctrl, delta, max_maf))                                   
        
