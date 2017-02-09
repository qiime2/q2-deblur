# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# test base, test data, and support code derived from q2-deblur

import unittest

import skbio
import biom
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt)
from q2_types.feature_data import DNAFASTAFormat

from q2_deblur import denoise_16S, denoise_other
from q2_deblur._denoise import _load_table, _hash_ids


# shamelessly adapted from q2-dada2
def _sort_seqs(seqs):
    return sorted(list(seqs), key=lambda x: x.metadata['id'])

# structure shamelessly adapted from q2-dada2
class TestDenoiseUtil(TestPluginBase):
    package = 'q2_deblur.tests'

    def setUp(self):
        super().setUp()
        self.table = self.get_data_path('expected/util')

    def test_load_table(self):
        exp = biom.example_table.copy()
        obs = _load_table(self.table)
        self.assertEqual(obs, exp)

    def test_hash_ids(self):
        table = _load_table(self.table)
        exp = {'O1': '00594a175ce5a58f286d91ca0a6f15a2',
               'O2': 'bfc60714f62f15e1e72e0b21454cc99e'}
        obs = _hash_ids(table)
        self.assertEqual(obs, exp)
        tab_obs = sorted(table.ids(axis='observation'))
        tab_exp = sorted(list(exp.values()))
        self.assertEqual(tab_obs, tab_exp)


class TestDenoise16S(TestPluginBase):
    package = 'q2_deblur.tests'

    def setUp(self):
        super().setUp()
        self.demux_seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path('sample_seqs_16S'), 'r')

    def test_defaults(self):
        exp_tab = biom.load_table(
            self.get_data_path('expected/16S-default.biom'))
        exp_rep_seqs = list(
            skbio.io.read(self.get_data_path('expected/16S-default.fasta'),
                          'fasta', constructor=skbio.DNA, lowercase='ignore'))
        for seq in exp_rep_seqs:
            del seq.metadata['description']

        obs_tab, rep_seqs = denoise_16S(self.demux_seqs, 100)

        rep_seqs = _sort_seqs(rep_seqs)
        exp_rep_seqs = _sort_seqs(exp_rep_seqs)

        self.assertEqual(obs_tab, exp_tab)
        self.assertEqual(rep_seqs, exp_rep_seqs)

    def test_all_reads_filtered(self):
        with self.assertRaisesRegex(ValueError, 'filter'):
            denoise_16S(self.demux_seqs, 10000)

    def test_bad_values_fail(self):
        # Just confirm that the machinery works, anything more specific is just
        # restating the _valid_inputs dict which is more declarative than a
        # unit-test anyways.
        with self.assertRaisesRegex(ValueError, 'trim_length'):
            denoise_16S(self.demux_seqs, -123)

        with self.assertRaisesRegex(ValueError, 'min_size'):
            denoise_16S(self.demux_seqs, 100, min_size=-1)


# structure shamelessly adapted from q2-dada2
class TestDenoiseOther(TestPluginBase):
    package = 'q2_deblur.tests'

    def setUp(self):
        super().setUp()
        self.demux_seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path('sample_seqs_other'), 'r')
        self.ref_ar = Artifact.load(
            self.get_data_path('../../assets/test_reference.qza'))
        self.ref = self.ref_ar.view(DNAFASTAFormat)

    def test_defaults(self):
        fp = self.get_data_path('expected/other-default.biom')
        exp_tab = biom.load_table(fp)
        exp_rep_seqs = list(
            skbio.io.read(self.get_data_path('expected/other-default.fasta'),
                          'fasta', constructor=skbio.DNA, lowercase='ignore'))
        for seq in exp_rep_seqs:
            del seq.metadata['description']

        obs_tab, rep_seqs = denoise_other(self.demux_seqs, self.ref, 100)

        rep_seqs = _sort_seqs(rep_seqs)
        exp_rep_seqs = _sort_seqs(exp_rep_seqs)

        self.assertEqual(rep_seqs, exp_rep_seqs)
        self.assertEqual(obs_tab, exp_tab)

    def test_all_reads_filtered(self):
        with self.assertRaisesRegex(ValueError, 'filter'):
            denoise_other(self.demux_seqs, self.ref, 10000)

    def test_bad_values_fail(self):
        # Just confirm that the machinery works, anything more specific is just
        # restating the _valid_inputs dict which is more declarative than a
        # unit-test anyways.
        with self.assertRaisesRegex(ValueError, 'trim_length'):
            denoise_other(self.demux_seqs, self.ref, -123)

        with self.assertRaisesRegex(ValueError, 'min_reads'):
            denoise_other(self.demux_seqs, self.ref, 100, min_reads=-1)


if __name__ == '__main__':
    unittest.main()
