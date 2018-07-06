#!/usr/bin/env python
"""
calc_and_compare_power.py

Calculate power of collapsing and adaptive collapsing tests from pval files
"""

import argparse
import os
import numpy as np
 
def run(power_file, alpha):
    cwd = os.getcwd()
    adaptive_list = []
    collapsing_list = []
    for f in os.listdir(cwd):
        if f.endswith('adaptive_pvals.computed.txt'):
            with open(f, 'r') as infile:
                adaptive_list.append(float(infile.readline().strip()))
        elif f.endswith('collapsing_pvals.computed.txt'):
            with open(f, 'r') as infile:
                collapsing_list.append(float(infile.readline().strip()))
    adaptive_list = np.array(adaptive_list)
    collapsing_list = np.array(collapsing_list)
    ad_reject = np.where(adaptive_list < alpha)[0].shape[0]
    col_reject = np.where(collapsing_list < alpha)[0].shape[0]
    ad_power = ad_reject / float(len(adaptive_list)) 
    col_power = col_reject / float(len(collapsing_list))
    with open(power_file, 'w') as outfile:
        outfile.write('adaptive_power\t{}\n'.format(ad_power))
        outfile.write('collapsing_power\t{}\n'.format(col_power))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--power_file', help='Specify name of output power file')
    parser.add_argument('--alpha', type=float, help='Specify alpha level')
    args = parser.parse_args()
    run(args.power_file, args.alpha)
