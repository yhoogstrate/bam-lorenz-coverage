#!/usr/bin/env python

import pysam
import subprocess
from tqdm import tqdm
import matplotlib.pyplot as plt
import tempfile
import os
from multiprocessing import Process
import time



class bamlorenzcoverage:
    def __init__(self):
        pass
  
    def bam_file_to_idx(self, bam_file):
        """
        Coverage plot needs the zero-statistic - i.e. the number of genomic bases not covered by reads
        """
        
        # nicer way to ctrl killing the child process first and not have hangs with ctrl c
        #https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python
        
        idx_observed = {}
        
        # FIFO stream / named pipe instead of actual file - saves humongous amounts of disk space for temp files
        tmp_filename = os.path.join(tempfile.mkdtemp() + '.fifo')
        os.mkfifo(tmp_filename)

        # I tried this with the Threading class but this often didnt parallelize
        parallel_thread = Process(target=pysam.samtools.depth, args=['-a',bam_file], kwargs={'save_stdout': tmp_filename})
        parallel_thread.start()

        fh = os.open(tmp_filename, os.O_RDONLY)
        
        chunk = os.read(fh, 1024*4).decode('ascii')
        while chunk:
            split = chunk.split('\n')

            if chunk[-1] == '\n':
                lines = split
                chunk = os.read(fh, 1024*4).decode('ascii')
            else:
                lines = split[:-1]
                chunk = split[-1] + os.read(fh, 1024*4).decode('ascii')

            for line in lines:
                if line: # last line is empty line ('')
                    depth = line.split('\t',2)[-1]

                    if not depth in idx_observed:
                        idx_observed[depth] = 1
                    else:
                        idx_observed[depth]  += 1

            

        idx_observed = {int(key): value for (key, value) in idx_observed.items()}

        os.remove(tmp_filename)
        return idx_observed

    def bam_file_to_idx_slow_and_mem_unsafe(self, bam_file):
        """
        Coverage plot needs the zero-statistic - i.e. the number of genomic bases not covered by reads
        """
        idx_observed = {}
        depth = ''
        status = 0
        for char in pysam.samtools.depth('-a', bam_file):
            if char == '\n':
                if not depth in idx_observed:
                    idx_observed[depth] = 1
                else:
                    idx_observed[depth]  += 1

                status = 0
                depth = ''
            else:
                if status == 2:
                    depth += char
                elif char == '\t':
                    status +=1 

        idx_observed = {int(key): value for (key, value) in idx_observed.items()}
        return idx_observed


    def bam_file_to_idx_mem_unsafe(self, bam_file):
        """
        Coverage plot needs the zero-statistic - i.e. the number of genomic bases not covered by reads
        """
        idx_observed = {}
        for line in tqdm(pysam.samtools.depth('-a', bam_file, split_lines=True)):
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

    def estimate_cumulative_coverage_curves(self, idx_observed):
        """
        IN:
        0: 100
        1: 70
        2: 15
        3: 50
        5: 25
        
        OUT:
        0: 100+70+15+50+25
        1: 70+15+50+25
        2: 15+50+25
        3: 50+25
        5: 25
        
        -----
        
        work top-down to save accumulated positions
        """
        cumulative_coverage_curves = {'minimum_coverage_depth':[], 'percentage_genome_covered':[]}

        idx_observed_cumulative = {}
        accumulation = 0
        for key in sorted(idx_observed, reverse = True):
            accumulation += idx_observed[key]
            idx_observed_cumulative[key] = accumulation
        
        idx_observed_cumulative_proportional = {}
        for key in idx_observed_cumulative:
            idx_observed_cumulative_proportional[key] = 100.0 * idx_observed_cumulative[key] / accumulation

        for min_depth in sorted(idx_observed_cumulative_proportional):
            cumulative_coverage_curves['minimum_coverage_depth'].append(min_depth)
            cumulative_coverage_curves['percentage_genome_covered'].append(idx_observed_cumulative_proportional[min_depth])
            
        return cumulative_coverage_curves

    def export_cumulative_coverage_curves(self, cumulative_coverage_curves, output_stream):
        output_stream.write("X_minimum_coverage_depth\tY_percentage_genome_covered\n")
        for n in range(len(cumulative_coverage_curves['minimum_coverage_depth'])):
            output_stream.write(str(cumulative_coverage_curves['minimum_coverage_depth'][n]) + "\t" + str(cumulative_coverage_curves['percentage_genome_covered'][n]) + "\n")

    def export_cumulative_coverage_plot(self, cumulative_coverage_curves, output_file, min_percentage_covered = 0.5):
        n = len([_ for _ in cumulative_coverage_curves['percentage_genome_covered'] if _ >=  min_percentage_covered])
        
        plt.plot(cumulative_coverage_curves['minimum_coverage_depth'][1:n], cumulative_coverage_curves['percentage_genome_covered'][1:n],'-bo')
        plt.xlabel('Minimum coverage depth')
        plt.ylabel('Percenatge genome covered (>= '+str(round(min_percentage_covered,1))+'%)')
        plt.savefig(output_file)
        plt.gcf().clear()

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
        cumu_frac_of_genome = 0
        cumulative_sequenced_bases = 0
        
        for depth in sorted(idx_observed, reverse = True):
            if depth != 0:
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
            if depth != 0:
                lorenz_curves['fraction_reads'].append(round(1.0 * lorenz_idx[depth][1] / total_sequenced_bases, 4))
                lorenz_curves['fraction_genome'].append(round(1.0 * lorenz_idx[depth][0] / total_covered_positions_of_genome,4))

        return lorenz_curves

    def export_lorenz_curves(self, lorenz_curves, output_stream):
        output_stream.write("X-fraction-sequenced-bases\tY-fraction-genome-covered\n")
        for n in range(len(lorenz_curves['fraction_reads'])):
            output_stream.write(str(lorenz_curves['fraction_reads'][n]) + "\t" + str(lorenz_curves['fraction_genome'][n]) + "\n")

    def export_lorenz_plot(self, lorenz_curves, output_file):
        plt.plot([0.0,1.0],[0.0,1.0],'k--')
        plt.plot(lorenz_curves['fraction_reads'], lorenz_curves['fraction_genome'],'-bo')
        plt.xlabel('Fraction sequenced bases')
        plt.ylabel('Fraction covered genome')
        plt.savefig(output_file)
        plt.gcf().clear()

