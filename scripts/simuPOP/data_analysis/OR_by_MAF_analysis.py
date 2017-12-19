#!/usr/bin/env python
"""
OR_analysis.py

For all .beta files, pull out the ORs to get the distribution
"""

import argparse
import os
 
def analyze(output_file, region_file):
    # make list of all (simuPOP) map files in directory 
    # (assuming these are replicates with same parameters)
    regions = []
    with open(region_file, 'r') as infile:
        for line in infile:
            line = line.strip().split()
            regions.append([int(line[0]),int(line[1])])
    sites = []
    for region in regions:
        sites.extend(range(region[0],region[1]+1))
    beta_list = []
    cwd = os.getcwd()
    for f in os.listdir(cwd):
        if f.endswith('beta'):
            beta_list.append(f)
    ORs = []
    MAFs = [] 
    for betaf in beta_list:
        with open(betaf, 'r') as infile:
            header = infile.readline()
            for line in infile:
                line = line.strip().split('\t')
                if int(line[0]) not in sites:
                    continue
                if float(line[3]) > 1.0:
                    ORs.append(float(line[3]))
                    MAFs.append(float(line[1]))
    with open(output_file, 'w') as outfile:
        outfile.write('OR\tMAF\treg\n')
        for OR, MAF in zip(ORs, MAFs):
            outfile.write('{}\t{}\t{}\n'.format(OR,MAF,region_file[:-4]))
        #outfile.write('\n'.join(map(str,ORs)))
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output_file', help='Specify name of output file')
    parser.add_argument('--region_file', help='Specify name of region file')
    args = parser.parse_args()
    analyze(args.output_file, args.region_file)
