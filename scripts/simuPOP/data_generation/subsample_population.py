#!/usr/bin/env python
"""
subsample_population.py

Draw a sample of individuals from a .aff file until desired case/control numbers
are reached
"""

import luigi
import os
import random
import sys
from time import sleep

class RequireFiles(luigi.ExternalTask):
    """ 
    Class for checking that necessary .mut files are present
    param aff: (string) name of .aff file to be sampled from
    """
    
    aff = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.aff)

class Sample(luigi.Task):
    """
    Class for sampling from a .aff file
    param aff: (string) name of .aff file to be sampled from
    param casen: (int) number of cases to sample
    param ctrln: (int) number of controls to sample
    """

    aff = luigi.Parameter()
    casen = luigi.IntParameter()
    ctrln = luigi.IntParameter()

    def requires(self):
        return RequireFiles(aff=self.aff)
    
    def output(self):
        aff_prefix = self.aff[:-3] 
        return luigi.LocalTarget('{prefix}smp'.format(prefix=aff_prefix))
    
    def run(self):
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
    param casen: (int) desired number of cases
    param ctrln: (int) desired number of ctrls
    """
    
    casen = luigi.IntParameter()
    ctrln = luigi.IntParameter()
    
    def requires(self):
        cwd = os.getcwd()
        aff_list = []
        for f in os.listdir(cwd):
            if f.endswith('.aff'):
                aff_list.append(f)
        for aff in aff_list:
            yield Sample(aff=aff, casen=self.casen, ctrln=self.ctrln)
