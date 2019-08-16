# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# test base, test data, and support code derived from q2-deblur

import unittest

import skbio
import biom
import pandas as pd
import pandas.util.testing as pdt
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt)
from q2_types.feature_data import DNAFASTAFormat

from q2_deblur import denoise_16S, denoise_other
from q2_deblur._format import STATS_HEADER
from q2_deblur._denoise import (_load_table, _hash_ids, _fasta_counts,
                                _read_fastq_seqs)


# shamelessly adapted from q2-dada2
def _sort_seqs(seqs):
    return sorted(list(seqs), key=lambda x: x.metadata['id'])


def _align_table(table, reference):
    t = table.sort_order(reference.ids())
    return t.sort_order(reference.ids(axis='observation'),
                        axis='observation')


# structure shamelessly adapted from q2-dada2
class TestDenoiseUtil(TestPluginBase):
    package = 'q2_deblur.tests'

    def setUp(self):
        super().setUp()
        self.table = self.get_data_path('expected/util')

    def test_fasta_counts(self):
        exp_unique, exp_count = 2, 123
        obs_unique, obs_count = _fasta_counts(self.get_data_path('./'),
                                              'test', 'fasta')
        self.assertEqual(obs_unique, exp_unique)
        self.assertEqual(obs_count, exp_count)

    def test_read_fastq_seqs(self):
        exp = 3
        obs = _read_fastq_seqs(self.get_data_path('test.fastq.gz'))
        self.assertEqual(obs, exp)

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

        obs_tab, rep_seqs, stats = denoise_16S(self.demux_seqs, 100)

        rep_seqs = _sort_seqs(rep_seqs)
        exp_rep_seqs = _sort_seqs(exp_rep_seqs)
        obs_tab = _align_table(obs_tab, exp_tab)

        self.assertEqual(obs_tab, exp_tab)
        self.assertEqual(rep_seqs, exp_rep_seqs)
        self.assertEqual(list(stats.columns), STATS_HEADER[1:])
        self.assertEqual(len(stats), 0)

    def test_left_trim_len(self):
        obs_tab, rep_seqs, stats = denoise_16S(self.demux_seqs, 110,
                                               left_trim_len=10)

        self.assertEqual(len(obs_tab.ids(axis='sample')), 16)
        self.assertEqual(len(obs_tab.ids(axis='observation')), 20)
        self.assertEqual(len(list(rep_seqs)), 20)
        self.assertEqual(len(stats.index), 0)

    def test_all_reads_filtered(self):
        with self.assertRaisesRegex(ValueError, 'filter.*10000'):
            denoise_16S(self.demux_seqs, 10000)

    def test_bad_values_fail(self):
        # Just confirm that the machinery works, anything more specific is just
        # restating the _valid_inputs dict which is more declarative than a
        # unit-test anyways.
        with self.assertRaisesRegex(ValueError, 'trim_length'):
            denoise_16S(self.demux_seqs, -123)

        with self.assertRaisesRegex(ValueError, 'min_size'):
            denoise_16S(self.demux_seqs, 100, min_size=-1)

    def test_with_stats(self):
        # manually assessed based on temp output
        #                            derep   dblr    art   chim  ref    miss
        exp_stats = [('L1S208', 100, 11, 69, 11, 64, 0, 0, 1, 2, 5, 46, 0, 0),
                     ('L1S257', 100, 12, 67, 12, 63, 0, 0, 0, 0, 4, 43, 0, 0),
                     ('L1S57',  100, 11, 60, 11, 58, 0, 0, 0, 0, 4, 39, 0, 0),
                     ('L1S76',  100, 11, 74, 10, 70, 0, 0, 0, 0, 3, 43, 0, 0),
                     ('L2S155', 100, 12, 40, 12, 40, 0, 0, 0, 0, 7, 29, 0, 0),
                     ('L2S175', 100, 11, 44, 11, 42, 0, 0, 0, 0, 6, 33, 0, 0),
                     ('L2S309', 100, 10, 38, 10, 38, 0, 0, 0, 0, 3, 23, 0, 0),
                     ('L2S357', 100,  8, 42,  8, 42, 0, 0, 0, 0, 4, 33, 0, 0),
                     ('L3S294', 100, 10, 33, 10, 33, 0, 0, 0, 0, 4, 18, 0, 0),
                     ('L3S313', 100, 12, 42, 12, 42, 0, 0, 0, 0, 5, 28, 0, 0),
                     ('L4S112', 100,  9, 36,  9, 36, 0, 0, 0, 0, 8, 34, 0, 0),
                     ('L4S63',  100,  9, 33,  9, 33, 0, 0, 0, 0, 3, 19, 0, 0),
                     ('L5S155', 100, 10, 44, 10, 44, 0, 0, 0, 0, 5, 32, 0, 0),
                     ('L5S174', 100, 13, 50, 13, 48, 0, 0, 0, 0, 4, 25, 0, 0),
                     ('L6S20',  100,  9, 45,  8, 43, 0, 0, 0, 0, 6, 39, 0, 0),
                     ('L6S68',  100, 14, 35, 14, 35, 0, 0, 0, 0, 5, 14, 0, 0)]

        exp_stats = pd.DataFrame(exp_stats, columns=STATS_HEADER)
        exp_stats.set_index('sample-id', inplace=True)

        _, _, obs_stats = denoise_16S(self.demux_seqs, 100, sample_stats=True)
        pdt.assert_frame_equal(obs_stats, exp_stats)

    def test_with_underscore_in_id(self):
        bad_seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path('sample_seqs_16S_underscore'), 'r')
        with self.assertRaisesRegex(
                ValueError, 'Deblur cannot.*L3S_313.'):
            denoise_16S(bad_seqs, 100)


class TestDenoiseOther(TestPluginBase):
    # structure shamelessly adapted from q2-dada2
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

        obs_tab, rep_seqs, stats = denoise_other(self.demux_seqs, self.ref,
                                                 100)

        rep_seqs = _sort_seqs(rep_seqs)
        exp_rep_seqs = _sort_seqs(exp_rep_seqs)
        obs_tab = _align_table(obs_tab, exp_tab)

        self.assertEqual(rep_seqs, exp_rep_seqs)
        self.assertEqual(obs_tab, exp_tab)
        self.assertEqual(list(stats.columns), STATS_HEADER[1:])
        self.assertEqual(len(stats), 0)

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
