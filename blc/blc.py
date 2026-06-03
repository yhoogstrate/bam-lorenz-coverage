#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""[License: GNU General Public License v3 (GPLv3)]
"""

import pysam
import subprocess
import functools
import warnings
from collections import defaultdict
import matplotlib.pyplot as plt
import tempfile
import os
from multiprocessing import Process


def deprecated(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        warnings.warn(f"{func.__name__} is deprecated", DeprecationWarning, stacklevel=2)
        return func(*args, **kwargs)
    return wrapper


class BamLorenzCoverage:
    READ_BUFFER_SIZE = 256 * 1024

    def __init__(self):
        pass

    def bam_file_to_idx(self, bam_file, region=None, bed_regions=None):
        """
        Coverage plot needs the zero-statistic - i.e. the number of genomic bases not covered by reads
        """

        # nicer way to ctrl killing the child process first and not have hangs with ctrl c
        # https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python

        total_investigated_genomic_positions = 0
        idx_observed = defaultdict(int)

        # FIFO stream / named pipe instead of actual file - saves humongous amounts of disk space for temp files
        tmp_filename = os.path.join(tempfile.mkdtemp() + '.fifo')

        try:
            os.mkfifo(tmp_filename)

            cmd = ['-a', bam_file]
            if region:
                cmd = ['-r', region] + cmd
            elif bed_regions:
                cmd = ['-b', bed_regions] + cmd

            # I tried this with the Threading class but this often didnt parallelize
            parallel_thread = Process(target=pysam.samtools.depth, args=cmd, kwargs={'save_stdout': tmp_filename})
            parallel_thread.start()

            fh = os.open(tmp_filename, os.O_RDONLY)
            try:
                chunk = os.read(fh, self.READ_BUFFER_SIZE).decode('ascii')
                while chunk:
                    split = chunk.split('\n')

                    if chunk[-1] == '\n':
                        lines = split
                        chunk = os.read(fh, self.READ_BUFFER_SIZE).decode('ascii')
                    else:
                        lines = split[:-1]
                        chunk = split[-1] + os.read(fh, self.READ_BUFFER_SIZE).decode('ascii')

                    for line in lines:
                        if line:  # last line is empty line ('')
                            total_investigated_genomic_positions += 1
                            idx_observed[int(line.rpartition('\t')[2])] += 1
            finally:
                os.close(fh)
                parallel_thread.terminate()
                parallel_thread.join()
        finally:
            os.remove(tmp_filename)

        return (dict(idx_observed), total_investigated_genomic_positions)

    @deprecated
    def bam_file_to_idx_slow_and_mem_unsafe(self, bam_file):
        """
        Coverage plot needs the zero-statistic - i.e. the number of genomic bases not covered by reads
        """
        idx_observed = {}
        depth = ''
        status = 0
        for char in pysam.samtools.depth('-a', bam_file):
            if char == '\n':
                if depth not in idx_observed:
                    idx_observed[depth] = 1
                else:
                    idx_observed[depth] += 1

                status = 0
                depth = ''
            else:
                if status == 2:
                    depth += char
                elif char == '\t':
                    status += 1

        idx_observed = {int(key): value for (key, value) in idx_observed.items()}
        return idx_observed

    @deprecated
    def bam_file_to_idx_mem_unsafe(self, bam_file):
        """
        Coverage plot needs the zero-statistic - i.e. the number of genomic bases not covered by reads
        """
        from tqdm import tqdm
        idx_observed = {}
        for line in tqdm(pysam.samtools.depth('-a', bam_file, split_lines=True)):
            depth = line.split('\t', 2)[-1]

            if depth not in idx_observed:
                idx_observed[depth] = 1
            else:
                idx_observed[depth] += 1

        idx_observed = {int(key): value for (key, value) in idx_observed.items()}
        return idx_observed

    @deprecated
    def bam_file_to_idx_old(self, bam_file):
        size_investigated_region = 0
        idx_observed = {}

        # skip -a
        with subprocess.Popen(['samtools', 'depth', '-a', bam_file], stdout=subprocess.PIPE, universal_newlines=True) as popen:
            for line in iter(popen.stdout.readline, ""):
                depth = line.strip('\n').split('\t')[-1]

                if depth not in idx_observed:
                    idx_observed[depth] = 1
                else:
                    idx_observed[depth] += 1
                size_investigated_region += 1

        idx_observed = {int(key): value for (key, value) in idx_observed.items()}

        return (idx_observed, size_investigated_region)

    @deprecated
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

                if depth not in idx_observed:
                    idx_observed[depth] = 1
                else:
                    idx_observed[depth] += 1
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
        cumulative_coverage_curves = {'minimum_coverage_depth': [], 'percentage_genome_covered': []}

        idx_observed_cumulative = {}
        accumulation = 0
        for key in sorted(idx_observed, reverse=True):
            accumulation += idx_observed[key]
            idx_observed_cumulative[key] = accumulation

        for min_depth in sorted(idx_observed_cumulative):
            cumulative_coverage_curves['minimum_coverage_depth'].append(min_depth)
            cumulative_coverage_curves['percentage_genome_covered'].append(100.0 * idx_observed_cumulative[min_depth] / accumulation)

        return cumulative_coverage_curves

    def export_cumulative_coverage_curves(self, cumulative_coverage_curves, output_stream):
        output_stream.write("X_minimum_coverage_depth\tY_percentage_genome_covered\n")
        for depth, pct in zip(cumulative_coverage_curves['minimum_coverage_depth'], cumulative_coverage_curves['percentage_genome_covered']):
            output_stream.write(f"{depth}\t{pct}\n")

    def export_cumulative_coverage_plot(self, cumulative_coverage_curves, output_file, min_percentage_covered=0.5):
        n = sum(1 for x in cumulative_coverage_curves['percentage_genome_covered'] if x >= min_percentage_covered)

        plt.plot(cumulative_coverage_curves['minimum_coverage_depth'][1:n], cumulative_coverage_curves['percentage_genome_covered'][1:n], '-bo')
        plt.xlabel('Minimum coverage depth')
        plt.ylabel(f'Percenatge genome covered (>= {round(min_percentage_covered, 1)}%)')
        plt.savefig(output_file)
        plt.gcf().clear()

    def estimate_lorenz_curves(self, idx_observed):
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
        "Using 3/13 of the sequenced bases, 1/7 of the covered genome is covered"
        (denominator is covered positions only, i.e. depth > 0; ensures curve ends at (1,1))

        """
        lorenz_idx = {}

        previous_depth = None
        cumu_positions = 0
        cumulative_sequenced_bases = 0

        for depth in sorted(idx_observed, reverse=True):
            if depth != 0:
                positions = idx_observed[depth]
                sequenced_bases = idx_observed[depth] * depth

                cumu_positions = positions
                cumulative_sequenced_bases = sequenced_bases
                if previous_depth:
                    cumu_positions += lorenz_idx[previous_depth][0]
                    cumulative_sequenced_bases += lorenz_idx[previous_depth][1]
                lorenz_idx[depth] = [cumu_positions, cumulative_sequenced_bases]

                previous_depth = depth

        total_covered_positions_of_genome = cumu_positions
        total_sequenced_bases = cumulative_sequenced_bases

        previous = [0, 0]
        top = 0
        denom = total_sequenced_bases * total_covered_positions_of_genome * 2

        lorenz_curves = {'fraction_reads': [0.0], 'fraction_genome': [0.0]}
        for depth in sorted(lorenz_idx, reverse=True):
            if depth != 0:
                lorenz_curves['fraction_reads'].append(round(1.0 * lorenz_idx[depth][1] / total_sequenced_bases, 4))
                lorenz_curves['fraction_genome'].append(round(1.0 * lorenz_idx[depth][0] / total_covered_positions_of_genome, 4))

                top += (lorenz_idx[depth][1] - previous[1]) * (lorenz_idx[depth][0] - previous[0])
                top += (lorenz_idx[depth][1] - previous[1]) * (previous[0]) * 2

                previous = lorenz_idx[depth]

        roc = 1.0 * top / denom

        lorenz_curves['roc'] = roc
        lorenz_curves['total_sequenced_bases'] = total_sequenced_bases
        lorenz_curves['total_covered_positions_of_genome'] = total_covered_positions_of_genome
        return lorenz_curves

    def export_lorenz_curves(self, lorenz_curves, output_stream):
        output_stream.write("X-fraction-sequenced-bases\tY-fraction-genome-covered\n")
        for reads, genome in zip(lorenz_curves['fraction_reads'], lorenz_curves['fraction_genome']):
            output_stream.write(f"{reads}\t{genome}\n")

    def export_lorenz_plot(self, lorenz_curves, output_file, sign_digits=3):
        plt.plot([0.0, 1.0], [0.0, 1.0], 'k--')
        plt.plot(lorenz_curves['fraction_reads'], lorenz_curves['fraction_genome'], '-bo')
        plt.text(0.0, 0.95, f"ROC={lorenz_curves['roc']:.{sign_digits}f}", fontsize=14)
        plt.xlabel('Fraction sequenced bases')
        plt.ylabel('Fraction covered genome')
        plt.savefig(output_file)
        plt.gcf().clear()
