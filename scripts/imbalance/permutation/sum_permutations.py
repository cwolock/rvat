#!/usr/bin/env python
"""
Reads sum matrix and summary file from generate_null_sum_matrix.py and generates an
empirical permutation-based Q-Q plot.                              
"""

from functools import partial
from scipy import stats
from utilities import *
import argparse
import ctypes
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import random
import seaborn as sns

sns.set_style('darkgrid')
 
def permute(counter, statuses):                                         
    """                                                                 
    Permute affection statuses and look up p-vals for this configuration
    :param counter: (int) iteration counter for multiprocessing to keep track
    :param statuses: (list) randomized list of case/control statuses         
    """                                                                 
    # randomize the statuses                                            
    perm_statuses = np.array(random.sample(statuses, len(statuses)))    
    case_indices = np.where(perm_statuses == '2')[0]                    
    ctrl_indices = np.where(perm_statuses == '1')[0] 
    n_case = len(case_indices)
    n_ctrl = len(ctrl_indices)                   
    # make contingency table and calculate pvals for each gene          
    pvals = []
    for row in range(col_matrix.shape[0]):
        case_sums = col_matrix[row,case_indices]
        ctrl_sums = col_matrix[row,ctrl_indices]
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
        if str(pval) == 'nan':
            pval = 1.0     
        pvals.append(pval)
    pvals = sorted(pvals)                                               
    return pvals

def read_summary_file(summary_file):
    """
    Read in summary and add gene info to dictionary
    :param summary_file: (string) name of summary file
    """
    pvals = []
    #genes = []
    with open(summary_file, 'r') as infile:
        header = infile.readline()
        first_line = infile.readline().strip().split(',')
        pvals.append(float(first_line[14]))
        ncase = int(first_line[7]) + int(first_line[8])
        nctrl = int(first_line[10]) + int(first_line[11])
        for line in infile:
            line = line.strip().split(',')
            pvals.append(float(line[14]))
            #genes.append(line[0])
    return (pvals, ncase, nctrl)#, genes)

def read_matrix_file(matrix_file, ngenes):
    """
    Read in matrix file
    :param matrix_file: (string) name of matrix file                  
    :param ngenes: (int) number of genes (or general collapsing units)
    """
    with open(matrix_file, 'r') as infile:
        samps = infile.readline().strip().split('\t')
        nsamps = len(samps) - 1
        col_matrix = np.zeros((ngenes, nsamps))
        for i, line in enumerate(infile):
            line = line.strip().split('\t')
            col_matrix[i,:] = line[1:]
    return col_matrix

def plot_qq(qq_file, exp, obs, lower, upper):
    """
    Create QQ-plot of observed and expected FET p-values, with 2.5%ile and 97.5%ile
    bounds.                                                                    
    Lambda = slope of regression on chi-sq transformed obs and exp p-values, after 
         removal of p-values of 1 and genome-wide sig. p-values                    
    :param qq_file: (string) name of QQ-plot file to be created                    
    :param exp: (np array) expected p-values                                       
    :param obs: (np array) observed p-values                                       
    :param lower: (np array) 2.5%ile exp p-values                                  
    :param upper: (np array) 97.5%ile exp p-values                                 
    """ 
    # change p-values > 1 to 1
    exp[exp > 1] = 1.0
    obs[obs > 1] = 1.0
    
    # have been having issue with pvalues of 0 with ultra ultra rare variant in cases
    tr_exp = exp[(exp > 0) & (obs > 0)]
    tr_obs = obs[(exp > 0) & (obs > 0)]
    tr_lower = lower[(exp > 0) & (obs > 0)]
    tr_upper = upper[(exp > 0) & (obs > 0)]
    upper = tr_upper
    lower = tr_lower

    # calculate lambda by regression method
    # first remove p-values of 1 and < genome-wide significance
    gws = 0.05 / len(obs)
    reg_exp = exp[(obs > gws) & (exp > gws) & (exp < 1) & (obs < 1)]
    reg_obs = obs[(obs > gws) & (exp > gws) & (exp < 1) & (obs < 1)]
    reg_exp = stats.chi2.ppf(1-reg_exp, 1)
    reg_obs = stats.chi2.ppf(1-reg_obs, 1)
    # for lstsq, explanatory vars must be in column form
    reg_exp = reg_exp[:,np.newaxis]
    # this least squares regression forces intercept to be 0
    slope = np.linalg.lstsq(reg_exp, reg_obs)[0][0]
    lambda_factor = slope
    lambda_factor = np.around(lambda_factor, 5)
    
    # transform p-values to log10                                                      
    exp = -np.log10(tr_exp)                                 
    obs = -np.log10(tr_obs)                  
    lower = -np.log10(lower)                                               
    upper = -np.log10(upper)                                                     
                                                                                       
    # set limits of qq plot axes                                                       
    axisMax_obs = np.ceil(max(obs))                                              
    axisMax_exp = np.ceil(max(exp))                                              
                                                                                       
    # initialize qqplot with axes and labels                                           
    fig = plt.figure(figsize=(12,12))
    plt.xlim([0, axisMax_exp])                                                         
    plt.xlabel('Expected -log10(p)', fontsize=20)
    plt.ylim([0, axisMax_obs])                                                         
    plt.ylabel('Observed -log10(p)', fontsize=20)
    plt.title('QQ Plot: Observed vs. expected p-values. Lambda = {l}'.format(  
              l = lambda_factor), fontsize=20)

    # change size of axis tick labels
    plt.tick_params(axis='both', which='major', labelsize=12)

    # plot the points, exp on x axis, obs on y axis                                
    dataAx = fig.add_subplot(1,1,1)                                                
    dataAx.plot(exp, obs, 'r.', label='_nolegend_', markersize=12)                    
    
    # plot a diagonal line for comparison                                          
    lineAx = fig.add_subplot(1,1,1)                                                
    lineAx.plot([0,max(axisMax_obs,axisMax_exp)], [0,max(axisMax_obs,axisMax_exp)],
                 'b-', label='_nolegend_')                                         
    uppAx = fig.add_subplot(1,1,1)                                                 
    uppAx.plot(exp, upper, 'g-')                                          
    lowAx = fig.add_subplot(1,1,1)                                                 
    lowAx.plot(exp, lower, 'y-')                                       
    plt.legend(['2.5th percentile of expected p-values',  
                '97.5th percentile of expected p-values'],
                loc=2)
    plt.tight_layout()                                                             
    plt.savefig(qq_file)                                                      
                                   
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('summary_file', 
                        help='Specify name of summary file')
    parser.add_argument('matrix_file', 
                        help='Specify name of matrix file') 
    parser.add_argument('output', help='Specify prefix for output files')
    parser.add_argument('--nperms', 
                        type=partial(valid_numerical_argument, arg_name='nperms',
                                     min_value=0, max_value=10000, arg_type=int),
                        default=1000,
                        help='Specify number of permutations')
    args = parser.parse_args()
    
    # name the output file
    qq_file = args.output + '.qqplot.pdf'
    
    # read summary file to get dictionary of gene-level info 
    ## num of cases, controls, and maximum number of qualified samples in any gene
    observed_pvals, ncase, nctrl = read_summary_file(args.summary_file)
    
    # calculate number of genes
    ngenes = len(observed_pvals)
    
    # read matrix file to build collapsing matrix
    col_matrix = read_matrix_file(args.matrix_file, ngenes)
    
    # calculate number of samples
    nsamps = col_matrix.shape[1]
    
    # create shared memory version of collapsing matrix
    shared_array_base1 = mp.Array(ctypes.c_float, ngenes * nsamps)
    shared_col_matrix = np.ctypeslib.as_array(shared_array_base1.get_obj())
    shared_col_matrix = shared_col_matrix.reshape(ngenes, nsamps)

    # perform permutations to generate expected p-values
    status_l = ['2'] * ncase
    status_l.extend(['1'] * nctrl)
    pool = mp.Pool(processes=12)
    permute_results = pool.map(partial(permute, statuses=status_l), range(args.nperms))
    pool.close()
    pool.join()

    # populate array of permuted p-values
    perm_pvals = np.ones((ngenes, args.nperms)) 
    for i, pvals in enumerate(permute_results):         
        perm_pvals[:,i] = pvals

    # calculate 2.5%ile and 97.5%ile of permuted pvals  
    bottom_perc = np.percentile(perm_pvals, 2.5, axis=1)
    top_perc = np.percentile(perm_pvals, 97.5, axis=1)

    # make arrays of ordered observed and expected p-values
    obs_pvals = np.array(sorted(observed_pvals))
    exp_pvals = np.mean(perm_pvals, axis=1)

    # generate qq plot
    plot_qq(qq_file, exp_pvals, obs_pvals, bottom_perc, top_perc)
