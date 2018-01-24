#!/usr/bin/env python
"""
combine_regions.py

Combines .mut files based on regions specified in a set of regions file. 
Give a range of positions to take from each selection coefficient file 
"""

import itertools
import luigi
import os
import random
import sys

class RequireFiles(luigi.ExternalTask):
    """ 
    Class for checking that necessary .smp files are present
    param f: (string) name of file
    """
    
    f = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.f)

class Splice(luigi.Task):
    """
    Class for splicing regions from .mut files
    param regions_file: (str) name of regions file
    param mut_list: (list) list of .mut files in directory
    """
    
    regions_file = luigi.Parameter()
    mut_list = luigi.ListParameter()

    def requires(self):
        f_list = list(self.mut_list) + [self.regions_file]
        for f in f_list:
            return RequireFiles(f=f)
    
    def output(self):
        prefix = '.'.join(self.regions_file.split('.')[0:2])
        name = '.'.join([prefix,'spl'])
        return luigi.LocalTarget(name)
    
    def run(self):
        # make dict of sel coeffs and sites to come from that file
        region_dict = {}
        with open(self.regions_file, 'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                start, end = int(line[0]), int(line[1])
                sites = range(start, end + 1)
                region_dict[line[2]] = sites
        # dict to keep track of sites for each person
        indiv_dict = {} 
        for f in self.mut_list:
            with open(f, 'r') as infile:
                # find the sel coeff that matches the .mut file
                sites = []
                for k, v in region_dict.iteritems():
                    if k == f.split('.')[0]: 
                        # record the sites to be included from this region
                        sites = v
                # get the relevant sites
                for i, line in enumerate(infile):
                    line = line.strip().split()
                    indiv = int(line[0])
                    if indiv not in indiv_dict:
                        indiv_dict[indiv] = [[],[]]
                    indiv_dict[indiv][i%2].extend(
                        [int(pos) for pos in line[1:] if 
                            int(pos) in sites])
        # sort the individual haplotypes
        for k, v in indiv_dict.iteritems():
            v[0] = sorted(list(set(v[0])))
            v[0] = map(str, v[0])
            v[1] = sorted(list(set(v[1])))
            v[1] = map(str, v[1])
        # sort the people
        ordered = sorted(indiv_dict.items(), key=lambda x:x[0])
        # write the output
        with self.output().open('w') as outfile:
            for indiv in ordered:
                outfile.write(str(indiv[0]) + '\t' + '\t'.join(indiv[1][0]) + '\n' + 
                              str(indiv[0]) + '\t' + '\t'.join(indiv[1][1]) + '\n')
             
class Parallelize(luigi.WrapperTask):
    """
    Class for parallelizing the region splicing task
    """
    
    def requires(self):
        # make list of .mut files for splicing and .regions.txt files
        cwd = os.getcwd()
        mut_list = [] # list of .mut files present in directory
        regions_list = [] # list of .regions.txt files present in directory
        for f in os.listdir(cwd):
            if f.endswith('.mut'):
                mut_list.append(f)
            elif f.endswith('regions.txt'):
                regions_list.append(f) 
        for regions_file in regions_list:
            yield Splice(regions_file=regions_file, mut_list=mut_list)
