#!/usr/bin/env python
"""
Make .mut file tab-delimited
"""
import argparse

def reformat(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.strip().split()
            outfile.write('\t'.join(line) + '\n')

if __name__ == '__main__':                                                  
    parser = argparse.ArgumentParser(                                       
        description=__doc__,                                                
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)             
    parser.add_argument('--input_file', help='Specify name of input file')
    parser.add_argument('--output_file', help='Specify name of output file')
    args = parser.parse_args()                                              
    reformat(args.input_file, args.output_file)                             
