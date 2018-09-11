#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:


import unittest
import os
import sys
import pysam


def sam_to_sorted_bam(sam_file, fixed_bam, tmp_dir='/tmp'):
    basename, ext = os.path.splitext(os.path.basename(sam_file))
    bam_file = os.path.join(tmp_dir, basename + '.unsorted.bam')

    print(bam_file, file=sys.stderr)

    # Seems like the file must exist in advance
    fh = open(bam_file, 'w')
    fh.close()

    # not a joke, pysam 0.15.0 requires to define output twice in order to write it to disk and not print it to stdout as well
    pysam.view('-o', bam_file, '-b', sam_file, save_stdout=bam_file)
    pysam.sort('-o', fixed_bam, '-l', '1', bam_file)
    pysam.index(fixed_bam)
    os.remove(bam_file)

    return bam_file


def main():
    unittest.main()
