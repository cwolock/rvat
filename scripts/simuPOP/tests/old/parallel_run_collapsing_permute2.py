#!/usr/bin/env python
"""
parallel_run_adaptive_collapsing_permute.py

Collects *.smp, *.spl, *regions.txt files and runs adaptive collapsing in parallel 
"""

import luigi
import os
from collapsing_permute2 import *

class RequireFiles(luigi.ExternalTask):
    """ 
    Class for checking that necessary files are present
    param smp: (string) name of .smp file
    param spl: (string) name of .spl file
    param reg: (string) name of .regions.txt file
    """
    
    smp = luigi.Parameter()
    spl = luigi.Parameter()
    reg = luigi.Parameter()

    def output(self):
        return (luigi.LocalTarget(self.smp), 
                luigi.LocalTarget(self.spl), 
                luigi.LocalTarget(self.reg))

class RunCollapsing(luigi.Task):
    """
    Class for running adaptive collapsing
    param smp: (string) name of .smp file
    param spl: (string) name of .spl file
    param reg: (string) name of .regions.txt file
    param nperms: (int) number of permutations
    """

    smp = luigi.Parameter()
    spl = luigi.Parameter()
    reg = luigi.Parameter()
    nperms = luigi.IntParameter()
    
    def requires(self):
        return RequireFiles(smp=self.smp, spl=self.spl, reg=self.reg)
    
    def output(self):
        prefix = self.spl[:-3] 
        return luigi.LocalTarget('{prefix}collapsing_pvals.txt'.format(prefix=prefix))
    
    def run(self):
        results = collapse(self.smp, self.spl, self.reg, self.nperms)
        with self.output().open('w') as outfile:
            outfile.write('\n'.join(map(str, results)) + '\n')
             
class Parallelize(luigi.WrapperTask):
    """
    Class for parallelizing the RunCollapsing task
    param nperms: (int) number of permutations
    """
    
    nperms = luigi.IntParameter()
    
    def requires(self):
        cwd = os.getcwd()
        smp_list = []
        spl_list = []
        reg_list = []
        for f in os.listdir(cwd):
            if f.endswith('.smp'):
                smp_list.append(f)
            elif f.endswith('.spl'):
                spl_list.append(f)
            elif f.endswith('regions.txt'):
                reg_list.append(f)
        smp_list = sorted(smp_list)
        spl_list = sorted(spl_list)
        reg_list = sorted(reg_list)
        for smp, spl, reg in zip(smp_list, spl_list, reg_list):
            yield RunCollapsing(smp=smp, spl=spl, reg=reg, nperms=self.nperms)
