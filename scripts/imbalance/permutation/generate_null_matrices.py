#!/usr/bin/env python
"""
generate_null_sum_matrix.py

Test for case/control enrichment under the null (no causal variants)
with imbalanced cohort size using weighted sum test
"""

from utilities import *
import argparse
import numpy as np
import random
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('darkgrid')

def run(gene_size, n_case, n_ctrl, delta, max_maf, iters, output_prefix):
    coll_summ_file = output_prefix + '.coll.summary.csv'
    coll_matrix_file = output_prefix + '.coll.matrix.txt'
    w_sum_summ_file = output_prefix + '.w_sum.summary.csv'
    w_sum_matrix_file = output_prefix + '.w_sum.matrix.txt'
    sum_summ_file = output_prefix + '.sum.summary.csv'
    sum_matrix_file = output_prefix + '.sum.matrix.txt'
    log_file = output_prefix + '.log'
    write_log(log_file, gene_size, n_case, n_ctrl, delta, max_maf, iters)
    with open(sum_summ_file, 'w') as sum_summ_out, \
         open(sum_matrix_file, 'w') as sum_matrix_out, \
         open(coll_summ_file, 'w') as coll_summ_out, \
         open(coll_matrix_file, 'w') as coll_matrix_out, \
         open(w_sum_summ_file, 'w') as w_sum_summ_out, \
         open(w_sum_matrix_file, 'w') as w_sum_matrix_out:
        sum_summ_out.write('Gene Name,Total Variant,Total SNV,Total Indel,'
                         'Unique Variant,Unique SNV,Unique Indel,'
                         'Qualified Case,Unqualified Case,Qualified Case Freq,'
                         'Qualified Ctrl,Unqualified Ctrl,Qualified Ctrl Freq,'
                         'Enriched Direction,FET\n')
        w_sum_summ_out.write('Gene Name,Total Variant,Total SNV,Total Indel,'
                         'Unique Variant,Unique SNV,Unique Indel,'
                         'Qualified Case,Unqualified Case,Qualified Case Freq,'
                         'Qualified Ctrl,Unqualified Ctrl,Qualified Ctrl Freq,'
                         'Enriched Direction,FET\n')
        coll_summ_out.write('Gene Name,Total Variant,Total SNV,Total Indel,'
                         'Unique Variant,Unique SNV,Unique Indel,'
                         'Qualified Case,Unqualified Case,Qualified Case Freq,'
                         'Qualified Ctrl,Unqualified Ctrl,Qualified Ctrl Freq,'
                         'Enriched Direction\n')
        sum_matrix_out.write('gene/sample\t' + 
                         '\t'.join(map(str, range(n_case + n_ctrl))) + '\n')
        w_sum_matrix_out.write('gene/sample\t' + 
                         '\t'.join(map(str, range(n_case + n_ctrl))) + '\n')
        coll_matrix_out.write('gene/sample\t' + 
                         '\t'.join(map(str, range(n_case + n_ctrl))) + '\n')
        for i in range(iters):
            cases, ctrls, variants, afs = generate_null_gts(gene_size, n_case, 
                                                            n_ctrl, delta, max_maf) 
            # collapsing
            qcase = len(set(np.where(cases > 0)[0]))
            qctrl = len(set(np.where(ctrls > 0)[0]))
            coll_summ_out.write('{gene},{tv},{ts},{ti},{uv},{us},{ui},'
                '{qc},{uqc},{qcf},{qct},{uqct},{qctf},{en}\n'.format(        
                gene=i,tv='NA',ts='NA',ti='NA',uv='NA',us='NA',              
                ui='NA',qc=qcase,uqc=n_case-qcase,qcf=qcase/float(n_case),   
                qct=qctrl,uqct=n_ctrl-qctrl,qctf=qctrl/float(n_ctrl),        
                en='NA'))                                                    
            matrix_line = (['1']*qcase + ['0']*(n_case-qcase) +                          
                           ['1']*qctrl + ['0']*(n_ctrl-qctrl))                           
            matrix_line_randomized = random.sample(matrix_line, len(matrix_line))        
            coll_matrix_out.write('{}\t'.format(i) + '\t'.join(matrix_line_randomized) + 
                                  '\n')
            
            # sum
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
            if str(pval) == 'nan':
                pval = 1.0
            sum_summ_out.write('{gene},{tv},{ts},{ti},{uv},{us},{ui},'
                            '{qc},{uqc},{qcf},{qct},{uqct},{qctf},{en},{fet}\n'.format(
                            gene=i,tv='NA',ts='NA',ti='NA',uv='NA',us='NA',
                            ui='NA',qc=qcase,uqc=n_case-qcase,qcf=qcase/float(n_case),
                            qct=qctrl,uqct=n_ctrl-qctrl,qctf=qctrl/float(n_ctrl),
                            en='NA',fet=pval))
            sum_matrix_out.write('{}\t'.format(i) + '\t'.join(map(str, X_vector)) + '\n')
            
            # weighted sum
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
            if str(pval) == 'nan':
                pval = 1.0
            w_sum_summ_out.write('{gene},{tv},{ts},{ti},{uv},{us},{ui},'
                            '{qc},{uqc},{qcf},{qct},{uqct},{qctf},{en},{fet}\n'.format(
                            gene=i,tv='NA',ts='NA',ti='NA',uv='NA',us='NA',
                            ui='NA',qc=qcase,uqc=n_case-qcase,qcf=qcase/float(n_case),
                            qct=qctrl,uqct=n_ctrl-qctrl,qctf=qctrl/float(n_ctrl),
                            en='NA',fet=pval))
            w_sum_matrix_out.write('{}\t'.format(i) + '\t'.join(map(str, X_vector)) + '\n')
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
