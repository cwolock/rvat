#!/usr/bin/env python
"""
OR_analysis.py

For all .beta files, pull out the ORs to get the distribution
"""

import argparse
import os
 
def analyze(output_file):
    # make list of all (simuPOP) map files in directory 
    # (assuming these are replicates with same parameters)
    beta_list = []
    reg_list = []
    cwd = os.getcwd()
    for f in os.listdir(cwd):
        if f.endswith('beta'):
            beta_list.append(f)
        elif f.endswith('regions.txt'):
            reg_list.append(f)
    beta_list = sorted(beta_list)
    reg_list = sorted(reg_list)
    ORs = []
    MAFs = [] 
    locations = []
    for betaf, regf in zip(beta_list, reg_list):
        regions = {}
        with open(regf, 'r') as infile:
             for line in infile:
                line = line.strip().split()
                regions[line[2]] = range(int(line[0]), int(line[1]))
        with open(betaf, 'r') as infile:
            #header = infile.readline()
            for line in infile:
                line = line.strip().split('\t')
                #if int(line[0]) not in sites:
                #    continue
                if float(line[3]) > 1.0:
                    ORs.append(float(line[3]))
                    MAFs.append(float(line[1]))
                    location = None
                    for k, v in regions.iteritems():
                        if int(line[0]) in v:
                            location = k
                    locations.append(location) 
    with open(output_file, 'w') as outfile:
        outfile.write('OR\tMAF\treg\n')
        for OR, MAF, location in zip(ORs, MAFs, locations):
            outfile.write('{}\t{}\t{}\n'.format(OR,MAF,location))
        #outfile.write('\n'.join(map(str,ORs)))
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output_file', help='Specify name of output file')
    args = parser.parse_args()
    analyze(args.output_file)
