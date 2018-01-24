#!/usr/bin/env python
"""
combine_regions.py


"""

import luigi
import os
import random
import sys
from time import sleep

class RequireFiles(luigi.ExternalTask):
    """ 
    Class for checking that necessary .mut files are present
    param mut: (string) name of .mut file to be sampled from
    """
    
    aff = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.aff)

class Sample(luigi.Task):
    """
    Class for sampling from a .mut file
    param mut: (string) name of .mut file to be sampled from
    param i: (int) iteration counter
    param subsamp_size: (int) desired subsample size (number of individuals) 
    param pop_size: (int) original population size (number of individuals)
    """

    aff = luigi.Parameter()
    i = luigi.IntParameter()
    #subsamp_size = luigi.IntParameter()
    casen = luigi.IntParameter()
    ctrln = luigi.IntParameter()
    #pop_size = luigi.IntParameter()

    def requires(self):
        return RequireFiles(aff=self.aff)
    
    def output(self):
        aff_prefix = self.aff[:-3] 
        return luigi.LocalTarget('{prefix}{i}.smp'.format(prefix=aff_prefix,i=self.i))
    
    def run(self):
        #indivs = range(1,self.pop_size+1)
        cases = set()
        ctrls = set()
        with open(self.aff, 'r') as infile:
            for line in infile:
                line = line.strip().split()
                if line[1] == '0':
                    ctrls.add(line[0])
                else: 
                    cases.add(line[0])
        if self.casen > len(cases) or self.ctrln > len(ctrls): sys.exit(
                         'Desired sample size too large for population.')
        casesamp = set(random.sample(cases, self.casen))
        ctrlsamp = set(random.sample(ctrls, self.ctrln))
        cohort = sorted(map(int, list(casesamp.union(ctrlsamp))))
        with self.output().open('w') as outfile:
            for indiv in cohort:
                outfile.write('{ID}\t{stat}\n'.format(
                    ID=indiv,stat='1' if str(indiv) in casesamp else '0'))
             
class Parallelize(luigi.WrapperTask):
    """
    Class for parallelizing the subsampling task
    param subsampsize: (int) desired subsample size (number of individuals)
    param popsize: (int) original population size (number of individuals)
    param niters: (int) number of times to subsample from each original population
    """
    
    #subsampsize = luigi.IntParameter()
    casen = luigi.IntParameter()
    ctrln = luigi.IntParameter()
    #popsize = luigi.IntParameter()
    niters = luigi.IntParameter()
    
    def requires(self):
        cwd = os.getcwd()
        aff_list = []
        for f in os.listdir(cwd):
            if f.endswith('.aff'):
                aff_list.append(f)
        for i in range(self.niters):
            for aff in aff_list:
                yield Sample(aff=aff, i=i, casen=self.casen, ctrln=self.ctrln)
