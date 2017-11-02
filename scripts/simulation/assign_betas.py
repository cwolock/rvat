#!/usr/bin/env python
"""
assign_betas.py

Assigns beta to each locus in .map file
"""

import luigi
import numpy as np
import os

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
    param region_l: (list) list of lists of loci and proportion of total length
    param map_f: (str) name of .map file to assign betas to
    """

    region_d = luigi.DictParameter()
    site_d = luigi.DictParameter()
    map_f = luigi.Parameter()
    
    def requires(self):
        return RequireFiles(f=self.map_f)
    
    def output(self):
        prefix = self.map_f[:-3] 
        return luigi.LocalTarget('{prefix}beta'.format(prefix=prefix))
    
    def run(self):
        sites = self.site_d
        regions = self.region_d
        with open(self.map_f, 'r') as infile:
            for line in infile:
                line = line.strip().split('\t') 
                locus = int(line[0])
                sites[locus][1] = float(line[2])
                region = sites[locus][0]
                seg = 1 if int(line[1]) > 0 else 0
                regions[region][1] += seg
        for region, value in regions.iteritems():
            # divide regional PAF by number of variant sites to get per-variant PAF
            value[2] = value[0]/float(value[1])
        ordered_sites = sorted(sites.items(), key=lambda x: x[0])
        with self.output().open('w') as outfile:
            for site in ordered_sites:
                if site[1][1] != 0:
                    # abs(beta) = ln(1+(eta/2*MAF))
                    print site[0], site[1][0], site[1][1]
                    print regions[site[1][0]]
                    beta = np.log(1+(regions[site[1][0]][2]/float((2*site[1][1]))))
                    OR = np.exp(beta)
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
    """
    
    regions = luigi.Parameter()
    length = luigi.IntParameter()   
 
    def requires(self):
        # make list of regions [[loci], frac of total length of region, # seg sites]
        region_d = {}
        site_d = {}
        with open(self.regions, 'r') as infile:
            for i, line in enumerate(infile):
                line = line.strip().split('\t')
                start = int(line[0])
                end = int(line[1])
                bases = end - start + 1
                frac = bases/float(self.length)
                loci = range(start, end+1)
                region_d[i] = [frac, 0, 0]
                for locus in loci:
                    site_d[locus] = [i,0, 0]
        cwd = os.getcwd()
        map_list = []
        for f in os.listdir(cwd):
            if f.endswith('.map'):
                map_list.append(f)
        for map_f in map_list:
            yield Assign(map_f=map_f, region_d=region_d, site_d=site_d)
