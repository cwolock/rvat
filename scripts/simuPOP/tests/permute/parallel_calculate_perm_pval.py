#!/usr/bin/env python
"""
parallel_run_adaptive_collapsing_permute.py

Collects *.smp, *.spl, *regions.txt files and runs adaptive collapsing in parallel 
"""

import luigi
import os
import numpy as np

class RequireFiles(luigi.ExternalTask):
    """ 
    Class for checking that necessary files are present
    param pfile: (string) name of pvals.txt file
    """
    
    pfile = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.pfile)

class CalculateP(luigi.Task):
    """
    Class for running adaptive collapsing
    param pfile: (string) name of pvals.txt file
    param nperms: (int) number of permutations
    """

    pfile = luigi.Parameter()
    nperms = luigi.IntParameter()
    
    def requires(self):
        return RequireFiles(pfile=self.pfile)
    
    def output(self):
        prefix = self.pfile[:-3] 
        return luigi.LocalTarget('{prefix}computed.txt'.format(prefix=prefix))
     
    def run(self):
        results = [[] for x in range(self.nperms + 1)]
        with open(self.pfile, 'r') as infile:
            for i, line in enumerate(infile):
                results[i % (self.nperms + 1)].append(float(line.strip()))
        min_p_list = [min(result) for result in results]
        min_p_list = np.array(min_p_list)
        # less than or less than/equal to???
        #more_extreme = np.where(min_p_list[1:] < min_p_list[0])[0].shape[0]
        more_extreme = np.where(min_p_list[1:] <= min_p_list[0])[0].shape[0]
        perm_p = more_extreme / float(self.nperms)
        #variance = perm_p * (1-perm_p) / float(self.nperms)
        with self.output().open('w') as outfile:
            outfile.write(str(perm_p) + '\n')
             
class Parallelize(luigi.WrapperTask):
    """
    Class for parallelizing the CalculateP task
    param nperms: (int) number of permutations
    """
    
    nperms = luigi.IntParameter()
    
    def requires(self):
        cwd = os.getcwd()
        pfile_list = []
        for f in os.listdir(cwd):
            if f.endswith('pvals.txt'):
                pfile_list.append(f)
        for pfile in pfile_list:
            yield CalculateP(pfile=pfile, nperms=self.nperms)
