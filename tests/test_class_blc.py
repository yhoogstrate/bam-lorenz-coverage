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
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(idx, {0: 392, 1: 108})

    def test_002_estimate_idx_from_bam(self):
        #           x x x x x
        #           x x x x x
        #           x x x x x
        # - - - - - - - - - - - - - - -

        # pct covered [  0 read ] =  15 / 15  = 100.0%
        # pct covered [ 1+ read ] =   5 / 15  =  33.3%
        # pct covered [ 2+ read ] =   5 / 15  =  33.3%
        # pct covered [ 3+ read ] =   5 / 15  =  33.3%

        idx = {0: 10, 3: 5}

        b = bamlorenzcoverage()
        cc = b.estimate_cumulative_coverage_curves(idx)

        # import sys
        # print(cc, file=sys.stderr)
        self.assertDictEqual(cc, {'minimum_coverage_depth': [0, 3], 'percentage_genome_covered': [100.0, 100.0 * 1.0 / 3.0]})

    def test_003_estimate_idx_from_bam(self):
        #           x x x x x
        #           x x x x x x x x x x
        # x x x x x x x x x x x x x x x
        # - - - - - - - - - - - - - - - - - - - - - - - - -

        # pct covered [  0 read ] =  25 / 25  = 100.0%
        # pct covered [ 1+ read ] =  15 / 25  =  60.0%
        # pct covered [ 2+ read ] =  10 / 15  =  40.0%
        # pct covered [ 3+ read ] =   5 / 15  =  20.0%

        idx = {0: 10, 1: 5, 2: 5, 3: 5}

        b = bamlorenzcoverage()
        cc = b.estimate_cumulative_coverage_curves(idx)

        # import sys
        # print(cc, file=sys.stderr)

        self.assertDictEqual(cc, {'minimum_coverage_depth': [0, 1, 2, 3], 'percentage_genome_covered': [100.0, 60.0, 40.0, 20.0]})

    def test_004_test_splice_junction(self):
        test_id = 'blc_004'

        input_file_sam = TEST_DIR + "test_" + test_id + ".sam"
        input_file_bam = T_TEST_DIR + "test_" + test_id + ".bam"

        sam_to_sorted_bam(input_file_sam, input_file_bam)

        b = bamlorenzcoverage()
        idx = b.bam_file_to_idx(input_file_bam)

        # print(idx, file=sys.stderr)
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(idx, {0: 392, 1: 108})

    def test_005_deletion(self):
        test_id = 'blc_005'

        input_file_sam = TEST_DIR + "test_" + test_id + ".sam"
        input_file_bam = T_TEST_DIR + "test_" + test_id + ".bam"

        sam_to_sorted_bam(input_file_sam, input_file_bam)

        b = bamlorenzcoverage()
        idx = b.bam_file_to_idx(input_file_bam)

        # print(idx, file=sys.stderr)
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(idx, {0: 392, 1: 108})

    def test_006_insertion(self):
        test_id = 'blc_006'

        input_file_sam = TEST_DIR + "test_" + test_id + ".sam"
        input_file_bam = T_TEST_DIR + "test_" + test_id + ".bam"

        sam_to_sorted_bam(input_file_sam, input_file_bam)

        b = bamlorenzcoverage()
        idx = b.bam_file_to_idx(input_file_bam)

        # print(idx, file=sys.stderr)
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(idx, {0: 400, 1: 100})

    def test_007_stacking(self):
        test_id = 'blc_007'

        input_file_sam = TEST_DIR + "test_" + test_id + ".sam"
        input_file_bam = T_TEST_DIR + "test_" + test_id + ".bam"

        sam_to_sorted_bam(input_file_sam, input_file_bam)

        b = bamlorenzcoverage()
        idx = b.bam_file_to_idx(input_file_bam)

        # print(idx, file=sys.stderr)
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(idx, {0: 372, 1: 48, 2: 80})

    def test_008_lorenz_01(self):
        #     x x x x
        # - - - - - - - - - -

        idx = {0: 6, 1: 4}

        b = bamlorenzcoverage()
        lc = b.estimate_lorenz_curves(idx)

        # print(idx, file=sys.stderr)
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(lc, {'fraction_genome': [0.0, 1.0], 'fraction_reads': [0.0, 1.0], 'roc': 0.5})

    def test_009_lorenz_02(self):
        # everything covered, is at least covered even densely
        #     x x x x
        #     x x x x
        # - - - - - - - - - -

        idx = {0: 6, 2: 4}

        b = bamlorenzcoverage()
        lc = b.estimate_lorenz_curves(idx)

        # print(idx, file=sys.stderr)
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(lc, {'fraction_genome': [0.0, 1.0], 'fraction_reads': [0.0, 1.0], 'roc': 0.5})

    def test_010_lorenz_03(self):
        # everything covered, is at least covered even densely
        #         x x x x x
        #   x x x x x
        # - - - - - - - - - -

        idx = {0: 6, 1: 6, 2: 2}

        b = bamlorenzcoverage()
        lc = b.estimate_lorenz_curves(idx)

        # print(idx, file=sys.stderr)
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(lc, {'fraction_genome': [0.0, 1.0 * 2 / 8, 1.0], 'fraction_reads': [0.0, 1.0 * 4 / 10, 1.0], 'roc': 0.425})

    def test_011_lorenz_03(self):
        # everything covered, is at least covered even densely
        #             x x x x x
        #       x x x x x
        # - - - - - - - - - - - - - -

        test_id = 'blc_011'

        input_file_sam = TEST_DIR + "test_" + test_id + ".sam"
        input_file_bam = T_TEST_DIR + "test_" + test_id + ".bam"

        sam_to_sorted_bam(input_file_sam, input_file_bam)

        b = bamlorenzcoverage()
        idx = b.bam_file_to_idx(input_file_bam)

        # print(idx, file=sys.stderr)
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(idx, {0: 6, 1: 6, 2: 2})

        lc = b.estimate_lorenz_curves(idx)

        # print(idx, file=sys.stderr)
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(lc, {'fraction_genome': [0.0, 1.0 * 2 / 8, 1.0], 'fraction_reads': [0.0, 1.0 * 4 / 10, 1.0], 'roc': 0.425})


if __name__ == '__main__':
    main()
