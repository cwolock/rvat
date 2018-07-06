#!/usr/bin/env python
"""
calculate_power.py

Calculate power of test from pvals
"""

import argparse
 
def run(pval_file, power_file, alpha):
    with open(pval_file, 'r') as infile, open(power_file, 'w') as outfile:
        reject, fail_to = 0, 0
        for line in infile:
            line = line.strip()
            pval = float(line)
            if pval <= alpha:
                reject += 1
            else:
                fail_to += 1
        power = reject / float(reject + fail_to)
        outfile.write(str(power))
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--pval_file', help='Specify name of pval file')
    parser.add_argument('--power_file', help='Specify name of output power file')
    parser.add_argument('--alpha', type=float, help='Specify alpha level')
    args = parser.parse_args()
    run(args.pval_file, args.power_file, args.alpha)
