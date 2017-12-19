#!/usr/bin/env python
"""
assign_disease.py

Assigns disease status to individuals in .spl file  based on betas
from .beta file until desired cohort sizes reached
"""

import luigi
import numpy as np
import os
from time import sleep
from scipy import stats

class RequireFiles(luigi.ExternalTask):
    """ 
    Class for checking that necessary .mut files are present
    param f: (string) name of file
    """
    
    f = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.f)

class Assign(luigi.Task):
    """
    Class for generating .map from .spl
    param length: (int) number of sites
    param casen: (int) desired number of cases
    param ctrln: (int) desired number of ctrls
    param spl: (string) name of .spl file
    param beta: (string) name of .beta file 
    """

    length = luigi.IntParameter()
    casen = luigi.IntParameter()
    ctrln = luigi.IntParameter()
    spl = luigi.Parameter()
    beta = luigi.Parameter()

    def requires(self):
        return [RequireFiles(f=self.spl), RequireFiles(f=self.beta)]
    
    def output(self):
        prefix = self.spl[:-3] 
        return luigi.LocalTarget('{prefix}aff'.format(prefix=prefix))
    
    def run(self):
        betas = np.array([0.0]* self.length)
        # get the beta for each locus
        with open(self.beta, 'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                betas[int(line[0])-1] = float(line[2])
        # determine disease status for each individual in .spl file
        cases = 0
        ctrls = 0
        with open(self.spl, 'r') as infile, self.output().open('w') as outfile:
            for i, line in enumerate(infile):
                line = line.strip().split('\t')
                indiv1 = line[0]
                seg_sites = line[1:]
                line = infile.next()
                line = line.strip().split('\t')
                indiv2 = line[0]
                if indiv1 != indiv2:
                    print "no match"
                seg_sites.extend(line[1:])
                gts = np.array([0.0]*self.length)
                for site in seg_sites:
                    gts[int(site)-1] += 1
                X_beta = np.dot(gts, betas)
                t = -2.94 + X_beta
                prob = np.exp(t)/float(1+np.exp(t))
                affect = stats.bernoulli.rvs(prob, size=1)[0]
                if affect == 1 and cases == self.casen and ctrls != self.ctrln:
                    continue
                elif affect == 0 and ctrls == self.ctrln and cases != self.casen:
                    continue
                elif cases == self.casen and ctrls == self.ctrln:
                    break
                hom_ref = np.where(gts == 0)[0].shape[0]
                het = np.where(gts == 1)[0].shape[0]
                hom_alt = np.where(gts == 2)[0].shape[0]
                if affect == 1:
                    cases += 1
                else:
                    ctrls += 1
                outfile.write('{}\t{}\n'.format(indiv1, affect))

class Parallelize(luigi.WrapperTask):
    """
    Class for parallelizing the Assign task
    param length: (int) number of sites
    param casen: (int) desired number of cases
    param ctrln: (int) desired number of ctrls
    """
    
    length = luigi.IntParameter()
    casen = luigi.IntParameter()
    ctrln = luigi.IntParameter()
    
    def requires(self):
        print "hello"
        cwd = os.getcwd()
        spl_list = []
        beta_list = []
        for f in os.listdir(cwd):
            if f.endswith('.spl'):
                spl_list.append(f)
            elif f.endswith('.beta'):
                beta_list.append(f)
        spl_list = sorted(spl_list)
        beta_list = sorted(beta_list)
        for spl, beta in zip(spl_list, beta_list):
            yield Assign(spl=spl, beta=beta, length=self.length, 
                         casen=self.casen, ctrln=self.ctrln)
