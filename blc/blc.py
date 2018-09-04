#!/usr/bin/env python

import pysam
import subprocess
from tqdm import tqdm


class bamlorenzcoverage:
    def __init__(self):
        pass
    
    def bam_file_to_idx(self, bam_file):
        idx_observed = {}
        for line in tqdm(pysam.samtools.depth(bam_file, split_lines=True)):
            depth = line.split('\t',2)[-1]
            
            if not depth in idx_observed:
                idx_observed[depth] = 1
            else:
                idx_observed[depth]  += 1

        idx_observed = {int(key): value for (key, value) in idx_observed.items()}
        return idx_observed

    def bam_file_to_idx_old(self, bam_file):
        size_investigated_region = 0
        idx_observed = {}

        # skip -a
        with subprocess.Popen(['samtools', 'depth', bam_file], stdout=subprocess.PIPE, universal_newlines=True) as popen:
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

    def estimatea_lorenz_curves(self, idx_observed):
        """
        In:
        
        
        Transformation:
          =====
      =======
    =====   =====
-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
0 0 1 2 2 2 3 2 1 0 0 0 0 0 0 0

stats table:

       |        Cumu   |  Seq.    Cumu
Depth  | Freq    Freq  | Bases   Bases
--------------------------------------- 
  0    |   9   1+4+2+9 |  0*9      0     neglected
--------------------------------------- 
  1    |   2    1+4+2  |  1*2     0+2
  2    |   4     1+4   |  2*4    0+2+8
  3    |   1      1    |  1*3   0+2+8+3


cumu_freq = 1+4+2+9 = 16
cumu_bases = 0+2+8+3 = 13

        
        
        Out:
        "Using 3/13 of the sequenced bases, 1/16 of the genome is covered"
        
        """
        lorenz_idx = {}
        
        previous_depth = None
        for depth in sorted(idx_observed, reverse = True):
            frac_of_genome = idx_observed[depth]
            sequenced_bases = idx_observed[depth] * depth
            
            cumu_frac_of_genome = frac_of_genome
            cumulative_sequenced_bases = sequenced_bases
            if previous_depth:
                cumu_frac_of_genome += lorenz_idx[previous_depth][0]
                cumulative_sequenced_bases += lorenz_idx[previous_depth][1]
            lorenz_idx[depth] = [cumu_frac_of_genome, cumulative_sequenced_bases]

            previous_depth = depth

        total_covered_positions_of_genome = cumu_frac_of_genome
        total_sequenced_bases = cumulative_sequenced_bases

        lorenz_curves = {'fraction_reads':[], 'fraction_genome':[]}
        for depth in sorted(lorenz_idx, reverse=True):
            lorenz_curves['fraction_reads'].append(round(1.0 * lorenz_idx[depth][1] / total_sequenced_bases, 4))
            lorenz_curves['fraction_genome'].append(round(1.0 * lorenz_idx[depth][0] / total_covered_positions_of_genome,4))

        return lorenz_curves

    def export_lorenz_curves(self, lorenz_curves, output_stream):
        output_stream.write("X-fraction-sequenced-bases\tY-fraction-genome-covered\n")
        for n in range(len(lorenz_curves['fraction_reads'])):
            output_stream.write(str(lorenz_curves['fraction_reads'][n]) + "\t" + str(lorenz_curves['fraction_genome'][n]) + "\n")
