#!/usr/bin/env python
"""
generate_map_file.py

Reads a .mut file (actually .spl for naming convenience) and produces a file
with number of individuals and freq of alt allele at each position
"""

import luigi
import os
import random
from time import sleep

class RequireFiles(luigi.ExternalTask):
    """ 
    Class for checking that necessary .mut files are present
    param spl: (string) name of .spl file
    """
    
    spl = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.spl)

class GenMap(luigi.Task):
    """
    Class for generating .map from .spl
    param length: (int) number of sites
    param pop: (int) population size
    param spl: (string) name of .spl file 
    """

    length = luigi.IntParameter()
    pop = luigi.IntParameter()
    spl = luigi.Parameter()

    def requires(self):
        return RequireFiles(spl=self.spl)
    
    def output(self):
        prefix = self.spl[:-3] 
        return luigi.LocalTarget('{prefix}map'.format(prefix=prefix))
    
    def run(self):
        site_d = {n: [0,0] for n in range(1,self.length+1)}
        with open(self.spl, 'r') as infile:
             for line in infile:
                line = line.strip().split('\t')
                seg_sites = line[1:]
                for site in seg_sites:
                    site_d[int(site)][0] += 1
        with self.output().open('w') as outfile:
            for k, v in site_d.iteritems():
                v[1] = v[0] / float(2*self.pop)
                outfile.write('{site}\t{alt}\t{freq}\n'.format(
                    site=k,alt=v[0], freq=v[1]))
             
class Parallelize(luigi.WrapperTask):
    """
    Class for parallelizing the GenMap task
    param length: (int) number of sites
    param popsize: (int) number of individuals
    """
    
    length = luigi.IntParameter()
    popsize = luigi.IntParameter()
    
    def requires(self):
        cwd = os.getcwd()
        spl_list = []
        for f in os.listdir(cwd):
            if f.endswith('.spl'):
                spl_list.append(f)
        for spl in spl_list:
            yield GenMap(spl=spl, length=self.length, pop=self.popsize)
