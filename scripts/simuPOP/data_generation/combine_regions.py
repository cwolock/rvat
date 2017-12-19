#!/usr/bin/env python
"""
combine_regions.py

Combines .mut files based on regions specific in regions file. 
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
    Class for splicing from a .smp files
    param files: (list) contains lists of [file_name, s] for all files to be spliced
    param pop: (string) population for output file name
    param bounds: (dict) {s: [[start, end],[start, end]...]}
    param prefix: (string) prefix of regions file, used for naming output files
    """

    files = luigi.ListParameter()
    pop = luigi.Parameter()
    bounds = luigi.DictParameter()
    prefix = luigi.Parameter()

    def requires(self):
        for f in self.files:
            return RequireFiles(f=f[0])
    
    def output(self):
        name = '.'.join([self.prefix,self.pop,'spl'])
        return luigi.LocalTarget(name)
    
    def run(self):
        indiv_dict = {}
        for f in self.files:
            with open(f[0], 'r') as infile:
                sites = []
                for region in self.bounds[f[1]]:
                    start = region[0]
                    end = region[1]
                    sites.extend(range(start, end + 1))
                #start = self.bounds[f[1]][0]
                #end = self.bounds[f[1]][1] 
                for i, line in enumerate(infile):
                    line = line.strip().split()
                    indiv = int(line[0])
                    if indiv not in indiv_dict:
                        indiv_dict[indiv] = [[],[]]
                    indiv_dict[indiv][i%2].extend(
                        [int(pos) for pos in line[1:] if 
                            int(pos) in sites])
        for k, v in indiv_dict.iteritems():
            v[0] = sorted(list(set(v[0])))
            v[0] = map(str, v[0])
            v[1] = sorted(list(set(v[1])))
            v[1] = map(str, v[1])
        ordered = sorted(indiv_dict.items(), key=lambda x:x[0])
        with self.output().open('w') as outfile:
            for indiv in ordered:
                outfile.write(str(indiv[0]) + '\t' + '\t'.join(indiv[1][0]) + '\n' + 
                              str(indiv[0]) + '\t' + '\t'.join(indiv[1][1]) + '\n')
             
class Parallelize(luigi.WrapperTask):
    """
    Class for parallelizing the region splicing task
    param regions : (string) name of regions file determining regional boundaries 
    """
    
    regions = luigi.Parameter()
    
    def requires(self):
        # keep name of regions file for naming output files
        prefix = self.regions[:-4]
        # make list of boundaries for each region to be spliced
        bounds = {}
        with open(self.regions, 'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                # FIX THIS IF MULTIPLE REGIONS FOR SAME SEL COEFFICIENT
                #bounds[line[2]] = [int(line[0]), int(line[1])]
                if line[2] not in bounds:
                    bounds[line[2]] = []
                bounds[line[2]].append([int(line[0]), int(line[1])])
        # group .mut files for splicing
        cwd = os.getcwd()
        mut_dict = {} # dict of .smp files present in directory
        sels = set() # set of selection coefficients present in directory
        pops = set() # set of population indices present in directory
        #subsamps = set() # set of subsample indices present in dictionary
        for f in os.listdir(cwd):
            # imagining file naming system like {sel}.{pop}.mut
            if f.endswith('.mut'):
                f_spl = f.strip().split('.')
                if f_spl[0] not in bounds:
                    continue
                sels.add(f_spl[0])
                pops.add(f_spl[1])
                #subsamps.add(f_spl[2])
                # sel coeff and pop for each file
                mut_dict[f] = [f_spl[0],f_spl[1]]
        if len(bounds) != len(sels):
            sys.exit('Invalid regions file')
        if ((len(mut_dict) % len(sels) != 0) or 
            (len(mut_dict) % len(pops) != 0)):# or
            #(len(smp_dict) % len(subsamps) != 0)):
            sys.exit('Invalid number of .mut files present')
        #pop_smp_combs = list(itertools.product(pops, subsamps))
        pop_dict = {pop: [] for pop in pops}
        # for each pop, make of dict of the files (and sel coeff) to combine
        for k, v in mut_dict.iteritems():
            pop_dict[v[1]].append([k,v[0]])
        for k, v in pop_dict.iteritems():
            yield Splice(files=v,pop=k,bounds=bounds,prefix=prefix)
