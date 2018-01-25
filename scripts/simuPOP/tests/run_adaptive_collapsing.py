#!/usr/bin/env python
"""
run_adaptive_collapsing.py

Run adaptive collapsing analysis on all .smp, .spl, .regions.txt files in directory and
calculate power/type I error
"""
from collapsing import *
import argparse
import os
 
def run(output_file, alpha):
    smps = []
    spls = []
    regs = []
    cwd = os.getcwd()
    for f in os.listdir(cwd):
        if f.endswith('smp'):
            smps.append(f)
        elif f.endswith('spl'):
            spls.append(f)
        elif f.endswith('regions.txt'):
            regs.append(f)
    smps = sorted(smps)
    spls = sorted(spls)
    regs = sorted(regs)
    reject = 0
    fail = 0
    for smp, spl in zip(smps, spls, regs):
        results = collapse(smp, spl)
        print results
        """
        if pval <= alpha:
            reject += 1
        else:
            fail += 1
    power = reject / float((reject + fail))  
    with open(output_file, 'w') as outfile:
        outfile.write(str(power))
        """
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output_file', help='Specify name of output file')
    parser.add_argument('--alpha', type=float, help='Specify alpha level')
    args = parser.parse_args()
    run(args.output_file, args.alpha)
