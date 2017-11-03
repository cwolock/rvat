#!/usr/bin/env python
"""
OR_analysis.py

For all .beta files, bins and calculates proportions for odds ratios
"""

import argparse
import os

def update(d, odds_r):
    if 1 <= odds_r < 1.1:
        d[(1,1.1)][0] += 1
    elif 1.1 <= odds_r < 1.2:
        d[(1.1,1.2)][0] += 1
    elif 1.2 <= odds_r < 1.3:
        d[(1.2,1.3)][0] += 1
    elif 1.3 <= odds_r < 1.4:
        d[(1.3,1.4)][0] += 1
    elif 1.4 <= odds_r < 1.5:
        d[(1.4,1.5)][0] += 1
    elif 1.5 <= odds_r < 1.6:
        d[(1.5,1.6)][0] += 1
    elif 1.6 <= odds_r < 2:
        d[(1.6,2)][0] += 1
    elif 2 <= odds_r < 3:
        d[(2,3)][0] += 1
    elif 3 <= odds_r < 5:
        d[(3,5)][0] += 1
    elif 5 <= odds_r < 100:
        d[(5,100)][0] += 1
    return d
     
def calc_frac(d, tot):
    for k, v in d.iteritems():
        v[1] = v[0]/float(tot)
    return d
 
def analyze(output_file):
    # make list of all beta files in directory 
    # (assuming these are replicates with same parameters)
    beta_list = []
    cwd = os.getcwd()
    for f in os.listdir(cwd):
        if f.endswith('beta'):
            beta_list.append(f)
    bins = [(1,1.1),(1.1,1.2),(1.2,1.3),(1.3,1.4),(1.4,1.5),
            (1.5,1.6),(1.6,2),(2,3),(3,5),(5,100)]
    rares = 0
    comms = 0
    r_dist = {b: [0,0] for b in bins}
    c_dist = {b: [0,0] for b in bins}
    for f in beta_list:
        with open(f, 'r') as infile:
            for line in infile:
                line = line.strip().split()
                MAF = float(line[1])
                OR = float(line[3])
                if 0 < MAF <= 0.01:
                    rares += 1
                    r_dist = update(r_dist, OR)   
                elif MAF > 0.01: 
                    comms += 1
                    c_dist = update(c_dist, OR)
    r_dist = calc_frac(r_dist, rares)
    c_dist = calc_frac(c_dist, comms)
    with open(output_file, 'w') as outfile:
        outfile.write('bin\trcount\tccount\trfrac\tcfrac\n')
        for k, v in r_dist.iteritems():
             outfile.write('{BIN}\t{rc}\t{cc}\t{rf}\t{cf}\n'.format(
                BIN=k,rc=v[0],cc=c_dist[k][0],rf=v[1],cf=c_dist[k][1]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output_file', help='Specify name of output file')
    args = parser.parse_args()
    analyze(args.output_file)
