#!/usr/bin/env python
"""
make_random_region_files.py

Given a file detailing region sizes and selection coefficients, generate regions
files with random coordinates meeting the desired criteria 
"""

import itertools
import luigi
import os
import random
import sys

class Randomize(luigi.Task):
    """
    Class for randomly generating regions files
    param sizes: (dict) dict of region sizes from different sel coeffs
    param perm: (int) permutation counter
    param prefix: (string) file prefix
    param length: (int) length of chromosome
    """

    sizes = luigi.DictParameter()
    perm = luigi.IntParameter()
    prefix = luigi.Parameter()
    length = luigi.IntParameter()
    
    def output(self):
        name = '.'.join([self.prefix,str(self.perm),'regions.txt'])
        return luigi.LocalTarget(name)
    
    def run(self):
        with self.output().open('w') as outfile:
            done = set()
            for k, v in self.sizes.iteritems():
                # this only handles one size per sel, can change later
                size = v[0]
                start = random.sample(range(1, self.length - size + 2), 1)[0]
                end = start + size - 1
                done.add(range(start, end + 1))
                outfile.write('{st}\t{en}\t{sel}\n'.format(st=start,en=end,sel=k))
             
class Parallelize(luigi.WrapperTask):
    """
    Class for parallelizing the region splicing task
    param config: (string) name of config file
    param length: (int) length of chromosome
    param perms: (int) number of permutations
    """
    
    config = luigi.Parameter()
    length = luigi.IntParameter()
    perms = luigi.IntParameter()
    
    def requires(self):
        # keep name of regions file for naming output files
        prefix = self.config[:-4]
        # make list of boundaries for each region to be spliced
        sizes = {}
        with open(self.config, 'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                # this should work multiple regions, same sel coeff, but I haven't tested
                if line[1] not in sizes:
                    sizes[line[1]] = []
                sizes[line[1]].append(int(line[0]))
        # make sure desired regions do not add up to longer than the chromosome length
        for k, v in sizes.iteritems():
            if sum(v) > self.length:
                sys.exit('Desired regions too long')
        for perm in range(self.perms):
            yield Randomize(sizes=sizes, perm=perm, prefix=prefix, length=self.length)
