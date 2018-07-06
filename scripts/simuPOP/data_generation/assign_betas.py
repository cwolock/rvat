#!/usr/bin/env python
"""
assign_betas.py

Assigns beta to each locus in .map file
"""

import luigi
import numpy as np
import os
import random

class RequireFiles(luigi.ExternalTask):
    """ 
    Class for checking that necessary files are present
    param f: (string) name of file
    """
    
    f = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.f)

class Assign(luigi.Task):
    """
    Class for generating .map from .spl
    param reg_f: (str) name of regions file
    param map_f: (str) name of .map file to assign betas to
    param site_d: (dict) dict of site info
    param delta: (float) fraction of variants that are causal (beta != 0)
    param epsilon: (float) fraction of causal variants with negative beta (protective)
    """

    reg_f = luigi.Parameter()
    map_f = luigi.Parameter()
    delta = luigi.FloatParameter()
    epsilon = luigi.FloatParameter()
    PAF = luigi.FloatParameter()
    mode = luigi.Parameter()
    ratio = luigi.Parameter() 
    
    def requires(self):
        return RequireFiles(f=self.map_f)
    
    def output(self):
        prefix = self.map_f[:-3] 
        return luigi.LocalTarget('{prefix}beta'.format(prefix=prefix))
    
    def run(self):
        # region_d[region #] = [reg. PAF, [seg sites], per-var PAF, 
        #                       length of region, selection coefficient]
        # site_d[coordinate] = [region, MAF, causal indicator, protective indicator] 
        region_d = {}
        site_d = {}
        with open(self.reg_f, 'r') as infile:
            for i, line in enumerate(infile):
                line = line.strip().split('\t')
                start = int(line[0])
                end = int(line[1])
                name = line[2]
                bases = end - start + 1
                if self.mode == 'bases':
                    bases = end - start + 1
                elif self.mode == 'coeff':
                    s = float(name.strip('sel'))
                    if s != 0:
                        #bases = 10000 * s
                        bases = s
                    else: #bases = 1
                        bases = 0.0001
                #elif self.mode = 
                # list of sites in region
                loci = range(start, end+1)
                # initialize site and region dict
                region_d[i] = [0, [], 0, bases]
                for locus in loci:
                    site_d[locus] = [i, 0, 0, 0]
        # if not doing regional variation in PAF assignments
        if self.ratio == 'uniform':
            region_d = {0: [0, [], 0, np.sum([v[3] for k, v 
                                                   in region_d.iteritems()])]}
            for k, v in site_d.iteritems():
                v[0] = 0
        # sort the region dict (allows use of numpy array below)
        ordered = sorted(region_d.items(), key=lambda x:x[0])
        # get lengths of all the regions
        segments = np.array([x[1][3] for x in ordered])
        # divide up the total PAF using one of several possible schemas
        # 1) PAF ratio of two regions is equal to ratio of sqrts of lengths
        if self.ratio == 'sqrt':
            factors = np.sqrt(segments)
        # 2) PAF ratio of two regions is equal to ratio of ln of lengths
        elif self.ratio == 'ln':
            factors = np.log(segments+1)
        elif self.ratio == 'expo':
            factors = 10**segments
        # 3) PAF ratio of two regions is equal to ratio of lengths
        elif self.ratio == 'linear' or self.ratio == 'uniform':
            factors = segments
        for i, factor in enumerate(factors):
            # calculate the inverse of the proportion of PAF that each region gets
            term = sum(factors/float(factor))
            # get regional PAF
            frac = self.PAF/float(term)
            reg = ordered[i][0]
            # add to dict
            region_d[reg][0] = frac
        sites = site_d
        regions = region_d
        with open(self.map_f, 'r') as infile:
            for line in infile:
                # each line is a region
                line = line.strip().split('\t') 
                locus = int(line[0])
                # record MAF of site
                sites[locus][1] = float(line[2])
                # increment the number of segregating sites in the region
                region = sites[locus][0]
                #seg = 1 if int(line[1]) > 0 else 0
                #regions[region][1] += seg
                if int(line[1]) > 0: 
                    regions[region][1].append(locus)
        for region, value in regions.iteritems():
            # if there are any variants in this regions
            if len(value[1]) != 0:
                # sample from list of seg sites to get causal sites
                #print np.ceil(len(value[1])*self.delta)
                causal = random.sample(value[1], 
                                       int(np.ceil(len(value[1])*self.delta)))
                # sample from list of causal sites to get protective sites
                protect = random.sample(causal, 
                                        int(np.ceil(len(causal)*self.epsilon))) 
                # assign sites to their respective groups
                for locus in causal:
                    sites[locus][2] = 1
                for locus in protect:
                    sites[locus][3] = 1
                # divide reg. PAF by # of CAUSAL variant sites to get per-variant PAF
                #value[2] = value[0]/float(value[1]) if float(value[1]) != 0 else 0
                value[2] = value[0]/float(len(causal)) if len(causal) > 0 else 0
        # sort sites by position
        ordered_sites = sorted(sites.items(), key=lambda x: x[0])
        with self.output().open('w') as outfile:
            for site in ordered_sites:
                # if site is causal
                #if site[1][1] != 0:
                if site[1][2] == 1:
                    # calculate beta
                    # abs(beta) = ln(1+(eta/2*MAF))
                    beta = np.log(1+(regions[site[1][0]][2]/float((2*site[1][1]))))
                    # if site is protective
                    if site[1][3] == 1:
                        beta = -beta
                    OR = np.exp(beta)
                # otherwise assign beta of 0
                else:
                    beta = 0
                    OR = 1
                outfile.write('{loc}\t{maf}\t{beta}\t{OR}\n'.format(
                        loc=site[0],maf=site[1][1],beta=beta,OR=OR))

class Parallelize(luigi.WrapperTask):
    """
    Class for parallelizing the Assign task
    param regions: (str) name of regions file determining regional boundaries
    param length: (int) length of region
    param PAF: (float) total PAF
    param delta: (float) fraction of variants that are causal
    param epsilon: (float) fraction of causal variants with negativ beta (protective)
    param mode: (str) use region length of seleciton coefficient to determine PAF dist
    param ratio: (str) how to split PAF between regions
    """
    
    #regions = luigi.Parameter()
    #length = luigi.IntParameter()
    PAF = luigi.FloatParameter()  
    mode = luigi.Parameter() 
    delta = luigi.FloatParameter()
    epsilon = luigi.FloatParameter()
    ratio = luigi.Parameter()
 
    def requires(self):
        with open('log.txt', 'w') as outfile:
            outfile.write('\n'.join(map(str, [self.PAF, self.mode, 
                                              self.delta, self.epsilon, self.ratio])))
        # make list of map files
        cwd = os.getcwd()
        map_list = []
        reg_list = []
        for f in os.listdir(cwd):
            if f.endswith('.map'):
                map_list.append(f)
            elif f.endswith('regions.txt'):
                reg_list.append(f)
        map_list = sorted(map_list)
        reg_list = sorted(reg_list)
        for map_f, reg_f in zip(map_list, reg_list):
            yield Assign(map_f=map_f, reg_f=reg_f, ratio=self.ratio, 
                delta=self.delta, epsilon=self.epsilon, PAF=self.PAF, mode=self.mode)
