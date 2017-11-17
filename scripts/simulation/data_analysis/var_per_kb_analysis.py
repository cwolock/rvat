#!/usr/bin/env python
"""
OR_analysis.py

For all .map files, calculates number of seg sites per kb
"""

import argparse
import os
 
def analyze(output_file, length):
    # make list of all (simuPOP) map files in directory 
    # (assuming these are replicates with same parameters)
    map_list = []
    cwd = os.getcwd()
    l = []
    for f in os.listdir(cwd):
        if f.endswith('map'):
            map_list.append(f)
    for mapf in map_list:
        with open(mapf, 'r') as infile:
            header = infile.readline()
            for i, line in enumerate(infile):
                pass
            varper = (1000*i)/float(length)
            l.append(varper)
    with open(output_file, 'w') as outfile:
        outfile.write('vpk\n')
        outfile.write('\n'.join(map(str,l)))
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output_file', help='Specify name of output file')
    parser.add_argument('--length', type=int, help='length of simulated segment')
    args = parser.parse_args()
    analyze(args.output_file, args.length)
