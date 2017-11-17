#!/usr/bin/env python
"""
combine_regions.py

Combines .smp files based on regions specific in regions file. 
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
    param comb: (tuple) (population, subsampl) tuple for output file name
    param bounds: (dict) {s: [start, end]}
    param prefix: (string) prefix of regions file, used for naming output files
    """

    files = luigi.ListParameter()
    comb = luigi.TupleParameter()
    bounds = luigi.DictParameter()
    prefix = luigi.Parameter()

    def requires(self):
        for f in self.files:
            return RequireFiles(f=f[0])
    
    def output(self):
        name = '.'.join([self.prefix,self.comb[0],self.comb[1],'spl'])
        return luigi.LocalTarget(name)
    
    def run(self):
        indiv_dict = {}
        for f in self.files:
            with open(f[0], 'r') as infile:
                start = self.bounds[f[1]][0]
                end = self.bounds[f[1]][1] 
                for i, line in enumerate(infile):
                    line = line.strip().split('\t')
                    indiv = int(line[0])
                    if indiv not in indiv_dict:
                        indiv_dict[indiv] = [[],[]]
                    indiv_dict[indiv][i%2].extend(
                        [int(pos) for pos in line[1:] if 
                            int(pos) >= start and int(pos) <= end])
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
    mode = luigi.Parameter()
    length = luigi.IntParameter()
    
    def requires(self):
        # keep name of regions file for naming output files
        prefix = self.regions[:-4]
        # make list of boundaries for each region to be spliced
        bounds = {}
        with open(self.regions, 'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                start = int(line[0])
                end = int(line[1])
                if self.mode == 'normal':
                    bounds[line[2]] = [int(line[0]), int(line[1])]
                elif self.mode == 'random':
                    dist = end - start + 1
                    new_start = random.randint(1,self.length-dist+2)
                    bounds[line[2]] = [new_start, new_start + dist - 1]
                    print new_start, new_start + dist - 1
        
        # group .muts files for splicing
        cwd = os.getcwd()
        smp_dict = {} # dict of .smp files present in directory
        sels = set() # set of selection coefficients present in directory
        pops = set() # set of population indices present in directory
        subsamps = set() # set of subsample indices present in dictionary
        for f in os.listdir(cwd):
            # imagining file naming system like {sel}.{pop}.{subsamp}.smp
            #if f.endswith('.smp'):
            if f.endswith('.mut'):
                f_spl = f.strip().split('.')
                if f_spl[0] not in bounds:
                    continue
                sels.add(f_spl[0])
                pops.add(f_spl[1])
                subsamps.add(f_spl[2])
                smp_dict[f] = [f_spl[0],f_spl[1], f_spl[2]]
        if len(bounds) != len(sels):
            sys.exit('Invalid regions file')
        if ((len(smp_dict) % len(sels) != 0) or 
            (len(smp_dict) % len(pops) != 0) or
            (len(smp_dict) % len(subsamps) != 0)):
            sys.exit('Invalid number of .smp files present')
        pop_smp_combs = list(itertools.product(pops, subsamps))
        comb_dict = {comb: [] for comb in pop_smp_combs}
        #comb_dict = {comb: [] for comb in 
        for k, v in smp_dict.iteritems():
            comb_dict[(v[1],v[2])].append([k,v[0]])
        for k, v in comb_dict.iteritems():
            yield Splice(files=v,comb=k,bounds=bounds,prefix=prefix)
