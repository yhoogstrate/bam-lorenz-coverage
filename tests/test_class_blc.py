#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:


import unittest
import os
from blc.blc import bamlorenzcoverage
from utils import main, sam_to_sorted_bam


TEST_DIR = "tests/blc/"
T_TEST_DIR = "tmp/" + TEST_DIR


# Nosetests doesn't use main()
if not os.path.exists(T_TEST_DIR):
    os.makedirs(T_TEST_DIR)


class TestIntronicBreakDetection(unittest.TestCase):
    def test_001_estimate_idx_from_bam(self):
        test_id = 'blc_001'

        input_file_sam = TEST_DIR + "test_" + test_id + ".sam"
        input_file_bam = T_TEST_DIR + "test_" + test_id + ".bam"

        sam_to_sorted_bam(input_file_sam, input_file_bam)

        b = bamlorenzcoverage()
        idx = b.bam_file_to_idx(input_file_bam)

        # print(idx, file=sys.stderr)
        self.assertDictEqual(idx, {0: 392, 1: 108})


if __name__ == '__main__':
    main()
