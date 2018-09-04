#!/usr/bin/env python


import subprocess
from tqdm import tqdm


class blc:
    def __init__(self):
        pass
    
    def bam_file_to_idx(self, bam_file):
        for line in tqdm(pysam.samtools.depth('-a', bam, split_lines=True)):
            depth = line.split('\t',2)[-1]
            
            if not depth in idx_observed:
                idx_observed[depth] = 1
            else:
                idx_observed[depth]  += 1
            size_investigated_region += 1

    def bam_file_to_idx_old(self, bam_file):
        size_investigated_region = 0
        idx_observed = {}

        with subprocess.Popen(['samtools', 'depth', '-a', bam_file], stdout=subprocess.PIPE, universal_newlines=True) as popen:
            for line in iter(popen.stdout.readline, ""):
                depth = line.strip('\n').split('\t')[-1]
                
                if not depth in idx_observed:
                    idx_observed[depth] = 1
                else:
                    idx_observed[depth]  += 1
                size_investigated_region += 1
                
        idx_observed = {int(key): value for (key, value) in idx_observed.items()}

        return (idx_observed, size_investigated_region)

    def coverage_file_to_idx(self, coverage_file):
        """
        do this initially from a coverage file - later pysam.depth direct parsing
        
        idx_observed = {depth: frequency}
        """

        size_investigated_region = 0
        idx_observed = {}

        with open(coverage_file, 'r') as fh:
            for line in fh:
                depth = line.strip('\n').split('\t')[-1]
                
                if not depth in idx_observed:
                    idx_observed[depth] = 1
                else:
                    idx_observed[depth]  += 1
                size_investigated_region += 1
                
        idx_observed = {int(key): value for (key, value) in idx_observed.items()}

        return (idx_observed, size_investigated_region)


    def estimate_cumulative_coverage_plot(self, idx_observed):
        """
        0: 100
        1: 70
        2: 15
        3: 50
        5: 25
        
        ->
        
        0: 100+70+15+50+25
        1: 70+15+50+25
        2: 15+50+25
        3: 50+25
        5: 25
        
        -----
        
        work top-down to save accumlated positions
        """
        idx_observed_cumulative = {}
        
        accumulation = 0
        for key in sorted(idx_observed, reverse = True):
            accumulation += idx_observed[key]
            idx_observed_cumulative[key] = accumulation
        
        idx_observed_cumulative_proportional = {}
        for key in idx_observed_cumulative:
            idx_observed_cumulative_proportional[key] = 100.0 * idx_observed_cumulative[key] / accumulation

        for min_depth in sorted(idx_observed_cumulative_proportional):
            print(str(min_depth) + "\t" + str(idx_observed_cumulative_proportional[min_depth]))

