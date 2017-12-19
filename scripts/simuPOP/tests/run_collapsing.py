#!/usr/bin/env python
"""
run_collapsing.py

Run collapsing analysis on all .smp and .spl files in directory and
calculate power/type I error
"""
from collapsing import *
import argparse
import os
 
def run(output_file, alpha):
    smps = []
    spls = []
    cwd = os.getcwd()
    for f in os.listdir(cwd):
        if f.endswith('smp'):
            smps.append(f)
        elif f.endswith('spl'):
            spls.append(f)
    smps = sorted(smps)
    spls = sorted(spls)
    reject = 0
    fail = 0
    #for smp, spl in zip(smps, spls):
    for smp in smps:
        # find .spl corresponding to .smp (multiple .smps go with one .spl)
        listed = smp.split('.')
        prefix1 = [listed[0], listed[1]]
        for f in spls:
            listed = f.split('.')
            prefix2 = [listed[0], listed[1]]
            if prefix1 == prefix2:
                spl = f
                break
        pval = collapse(smp, spl)
        if pval <= alpha:
            reject += 1
        else:
            fail += 1
    power = reject / float((reject + fail))  
    with open(output_file, 'w') as outfile:
        outfile.write(str(power))
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output_file', help='Specify name of output file')
    parser.add_argument('--alpha', type=float, help='Specify alpha level')
    args = parser.parse_args()
    run(args.output_file, args.alpha)
