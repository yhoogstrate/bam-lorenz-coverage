#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""[License: GNU General Public License v3 (GPLv3)]
"""


import unittest
import os
import io
import tempfile
import re
from blc.blc import BamLorenzCoverage
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

        b = BamLorenzCoverage()
        idx, n = b.bam_file_to_idx(input_file_bam)

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

        b = BamLorenzCoverage()
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

        b = BamLorenzCoverage()
        cc = b.estimate_cumulative_coverage_curves(idx)

        # import sys
        # print(cc, file=sys.stderr)

        self.assertDictEqual(cc, {'minimum_coverage_depth': [0, 1, 2, 3], 'percentage_genome_covered': [100.0, 60.0, 40.0, 20.0]})

    def test_004_test_splice_junction(self):
        test_id = 'blc_004'

        input_file_sam = TEST_DIR + "test_" + test_id + ".sam"
        input_file_bam = T_TEST_DIR + "test_" + test_id + ".bam"

        sam_to_sorted_bam(input_file_sam, input_file_bam)

        b = BamLorenzCoverage()
        idx, n = b.bam_file_to_idx(input_file_bam)

        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(idx, {0: 392, 1: 108})

        # additional stats
        self.assertEqual(n, 500)

    def test_005_deletion(self):
        test_id = 'blc_005'

        input_file_sam = TEST_DIR + "test_" + test_id + ".sam"
        input_file_bam = T_TEST_DIR + "test_" + test_id + ".bam"

        sam_to_sorted_bam(input_file_sam, input_file_bam)

        b = BamLorenzCoverage()
        idx, n = b.bam_file_to_idx(input_file_bam)

        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(idx, {0: 392, 1: 108})

    def test_006_insertion(self):
        test_id = 'blc_006'

        input_file_sam = TEST_DIR + "test_" + test_id + ".sam"
        input_file_bam = T_TEST_DIR + "test_" + test_id + ".bam"

        sam_to_sorted_bam(input_file_sam, input_file_bam)

        b = BamLorenzCoverage()
        idx, n = b.bam_file_to_idx(input_file_bam)

        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(idx, {0: 400, 1: 100})

        # additional stats
        self.assertEqual(n, 500)

    def test_007_stacking(self):
        test_id = 'blc_007'

        input_file_sam = TEST_DIR + "test_" + test_id + ".sam"
        input_file_bam = T_TEST_DIR + "test_" + test_id + ".bam"

        sam_to_sorted_bam(input_file_sam, input_file_bam)

        b = BamLorenzCoverage()
        idx, n = b.bam_file_to_idx(input_file_bam)

        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(idx, {0: 372, 1: 48, 2: 80})

        # additional stats
        self.assertEqual(n, 500)

    def test_008_lorenz_01(self):
        #     x x x x
        # - - - - - - - - - -

        idx = {0: 6, 1: 4}

        b = BamLorenzCoverage()
        lc = b.estimate_lorenz_curves(idx)

        # print(idx, file=sys.stderr)
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertListEqual(lc['fraction_genome'], [0.0, 1.0])
        self.assertListEqual(lc['fraction_reads'], [0.0, 1.0])

        # additional stats
        self.assertEqual(lc['roc'], 0.5)

    def test_009_lorenz_02(self):
        # everything covered, is at least covered even densely
        #     x x x x
        #     x x x x
        # - - - - - - - - - -

        idx = {0: 6, 2: 4}

        b = BamLorenzCoverage()
        lc = b.estimate_lorenz_curves(idx)

        # print(idx, file=sys.stderr)
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertListEqual(lc['fraction_genome'], [0.0, 1.0])
        self.assertListEqual(lc['fraction_reads'], [0.0, 1.0])

        # additional stats
        self.assertEqual(lc['roc'], 0.5)

    def test_010_lorenz_03(self):
        # everything covered, is at least covered even densely
        #         x x x x x
        #   x x x x x
        # - - - - - - - - - -

        idx = {0: 6, 1: 6, 2: 2}

        b = BamLorenzCoverage()
        lc = b.estimate_lorenz_curves(idx)

        # print(idx, file=sys.stderr)
        self.assertListEqual(lc['fraction_genome'], [0.0, 1.0 * 2 / 8, 1.0])
        self.assertListEqual(lc['fraction_reads'], [0.0, 1.0 * 4 / 10, 1.0])

        # additional stats
        self.assertEqual(lc['roc'], 0.425)

    def test_011_lorenz_03(self):
        # everything covered, is at least covered even densely
        #             x x x x x
        #       x x x x x
        # - - - - - - - - - - - - - -

        test_id = 'blc_011'

        input_file_sam = TEST_DIR + "test_" + test_id + ".sam"
        input_file_bam = T_TEST_DIR + "test_" + test_id + ".bam"

        sam_to_sorted_bam(input_file_sam, input_file_bam)

        b = BamLorenzCoverage()
        idx, n = b.bam_file_to_idx(input_file_bam)
        lc = b.estimate_lorenz_curves(idx)

        # print(idx, file=sys.stderr)
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertDictEqual(idx, {0: 6, 1: 6, 2: 2})

        # print(idx, file=sys.stderr)
        # denote that it only considers the (size of the) sequences described in the SAM header
        self.assertListEqual(lc['fraction_genome'], [0.0, 1.0 * 2 / 8, 1.0])
        self.assertListEqual(lc['fraction_reads'], [0.0, 1.0 * 4 / 10, 1.0])

        # additional stats
        self.assertEqual(n, 14)  # sam header say reference size is 14
        self.assertEqual(lc['roc'], 0.425)
        self.assertEqual(lc['total_sequenced_bases'], 10)
        self.assertEqual(lc['total_covered_positions_of_genome'], 8)

    def test_012_region(self):
        # everything covered, is at least covered even densely
        #             x x x x x
        #       x x x x x
        # - - - - - - - - - - - - - -

        test_id = 'blc_012'

        input_file_sam = TEST_DIR + "test_" + test_id + ".sam"
        input_file_bam = T_TEST_DIR + "test_" + test_id + ".bam"

        sam_to_sorted_bam(input_file_sam, input_file_bam)

        b = BamLorenzCoverage()
        idx, n = b.bam_file_to_idx(input_file_bam, 'chr1:2-14')

        self.assertEqual(n, 13)  # sam header say reference size is 14, but we start at 2nd position

    def test_013_bed(self):
        # everything covered, is at least covered even densely
        #             x x x x x
        #       x x x x x
        # - - - - - - - - - - - - - -

        test_id = 'blc_013'

        input_file_sam = TEST_DIR + "test_" + test_id + ".sam"
        input_file_bed = TEST_DIR + "test_" + test_id + ".bed"
        input_file_bam = T_TEST_DIR + "test_" + test_id + ".bam"

        sam_to_sorted_bam(input_file_sam, input_file_bam)

        b = BamLorenzCoverage()
        idx, n = b.bam_file_to_idx(input_file_bam, None, input_file_bed)

        self.assertEqual(n, 12)

    def test_014_export_cumulative_coverage_curves_tsv(self):
        idx = {0: 10, 1: 5, 2: 5}
        b = BamLorenzCoverage()
        cc = b.estimate_cumulative_coverage_curves(idx)

        output = io.StringIO()
        b.export_cumulative_coverage_curves(cc, output)

        lines = output.getvalue().strip().split('\n')
        self.assertEqual(lines[0], "X_minimum_coverage_depth\tY_percentage_genome_covered")
        self.assertEqual(len(lines), 4)
        self.assertTrue(all('\t' in line for line in lines[1:]))

    def test_015_export_cumulative_coverage_to_file(self):
        idx = {0: 10, 1: 5}
        b = BamLorenzCoverage()
        cc = b.estimate_cumulative_coverage_curves(idx)

        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.tsv') as f:
            temp_path = f.name

        try:
            with open(temp_path, 'w') as f:
                b.export_cumulative_coverage_curves(cc, f)

            with open(temp_path, 'r') as f:
                content = f.read()

            self.assertIn("X_minimum_coverage_depth", content)
            self.assertIn("Y_percentage_genome_covered", content)
        finally:
            os.unlink(temp_path)

    def test_016_export_lorenz_curves_tsv(self):
        idx = {0: 6, 1: 4}
        b = BamLorenzCoverage()
        lc = b.estimate_lorenz_curves(idx)

        output = io.StringIO()
        b.export_lorenz_curves(lc, output)

        lines = output.getvalue().strip().split('\n')
        self.assertEqual(lines[0], "X-fraction-sequenced-bases\tY-fraction-genome-covered")
        self.assertEqual(len(lines), 3)

    def test_017_empty_coverage_dict(self):
        idx = {}
        b = BamLorenzCoverage()

        try:
            cc = b.estimate_cumulative_coverage_curves(idx)
            self.assertEqual(len(cc['minimum_coverage_depth']), 0)
        except (ValueError, ZeroDivisionError):
            pass

    def test_018_single_depth_value(self):
        idx = {5: 10}
        b = BamLorenzCoverage()
        cc = b.estimate_cumulative_coverage_curves(idx)

        self.assertEqual(len(cc['minimum_coverage_depth']), 1)
        self.assertEqual(cc['percentage_genome_covered'][0], 100.0)

    def test_019_very_high_coverage(self):
        idx = {0: 1000, 1000: 100, 10000: 10}
        b = BamLorenzCoverage()
        cc = b.estimate_cumulative_coverage_curves(idx)

        self.assertEqual(cc['percentage_genome_covered'][0], 100.0)
        self.assertTrue(all(0 <= p <= 100 for p in cc['percentage_genome_covered']))

    def test_020_lorenz_with_single_depth(self):
        idx = {0: 90, 5: 10}
        b = BamLorenzCoverage()
        lc = b.estimate_lorenz_curves(idx)

        self.assertIsInstance(lc['roc'], float)
        self.assertGreaterEqual(lc['roc'], 0.0)
        self.assertLessEqual(lc['roc'], 1.0)

    def test_021_export_cumulative_coverage_plot_creates_file(self):
        idx = {0: 10, 1: 5, 2: 3}
        b = BamLorenzCoverage()
        cc = b.estimate_cumulative_coverage_curves(idx)

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.svg') as f:
            temp_path = f.name

        try:
            b.export_cumulative_coverage_plot(cc, temp_path)

            self.assertTrue(os.path.exists(temp_path))
            self.assertGreater(os.path.getsize(temp_path), 0)

            with open(temp_path, 'r') as f:
                content = f.read()
            self.assertIn('<svg', content)
            self.assertIn('</svg>', content)
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)

    def test_022_export_lorenz_plot_creates_file(self):
        idx = {0: 6, 1: 4}
        b = BamLorenzCoverage()
        lc = b.estimate_lorenz_curves(idx)

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.svg') as f:
            temp_path = f.name

        try:
            b.export_lorenz_plot(lc, temp_path)

            self.assertTrue(os.path.exists(temp_path))

            with open(temp_path, 'r') as f:
                content = f.read()

            self.assertIn('ROC=', content)
            self.assertIn('<svg', content)
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)

    def test_023_lorenz_plot_with_custom_precision(self):
        idx = {0: 6, 1: 4}
        b = BamLorenzCoverage()
        lc = b.estimate_lorenz_curves(idx)

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.svg') as f:
            temp_path = f.name

        try:
            b.export_lorenz_plot(lc, temp_path, sign_digits=1)

            with open(temp_path, 'r') as f:
                content = f.read()

            roc_match = re.search(r'ROC=(\d+\.\d+)', content)
            self.assertIsNotNone(roc_match)
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)


if __name__ == '__main__':
    main()
