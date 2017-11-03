#!/usr/bin/env python
"""
combine_regions.py


"""

import luigi
import os
import random
from time import sleep

class RequireFiles(luigi.ExternalTask):
    """ 
    Class for checking that necessary .mut files are present
    param mut: (string) name of .mut file to be sampled from
    """
    
    mut = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.mut)

class Sample(luigi.Task):
    """
    Class for sampling from a .mut file
    param mut: (string) name of .mut file to be sampled from
    param i: (int) iteration counter
    param subsamp_size: (int) desired subsample size (number of individuals) 
    param pop_size: (int) original population size (number of individuals)
    """

    mut = luigi.Parameter()
    i = luigi.IntParameter()
    subsamp_size = luigi.IntParameter()
    pop_size = luigi.IntParameter()

    def requires(self):
        return RequireFiles(mut=self.mut)
    
    def output(self):
        mut_prefix = self.mut[:-3] 
        return luigi.LocalTarget('{prefix}{i}.smp'.format(prefix=mut_prefix,i=self.i))
    
    def run(self):
        indivs = range(1,self.pop_size+1)
        subsamp = set(random.sample(indivs, self.subsamp_size))
        with open(self.mut) as infile, self.output().open('w') as outfile:
            n = 1
            for i, line in enumerate(infile):
                line = line.strip().split()
                if int(line[0]) in subsamp:
                    outfile.write(str(n) + '\t' + '\t'.join(line[1:]) + '\n')
                    if i%2 == 1:
                        n += 1
             
class Parallelize(luigi.WrapperTask):
    """
    Class for parallelizing the subsampling task
    param subsampsize: (int) desired subsample size (number of individuals)
    param popsize: (int) original population size (number of individuals)
    param niters: (int) number of times to subsample from each original population
    """
    
    subsampsize = luigi.IntParameter()
    popsize = luigi.IntParameter()
    niters = luigi.IntParameter()
    
    def requires(self):
        cwd = os.getcwd()
        mut_list = []
        for f in os.listdir(cwd):
            if f.endswith('.mut'):
                mut_list.append(f)
        for i in range(self.niters):
            for mut in mut_list:
                yield Sample(mut=mut, i=i, subsamp_size=self.subsampsize, 
                             pop_size=self.popsize)
