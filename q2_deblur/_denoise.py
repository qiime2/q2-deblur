# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import tempfile
import subprocess
import os
import hashlib

import biom
import skbio

from q2_types.per_sample_sequences import \
        SingleLanePerSampleSingleEndFastqDirFmt
from q2_types.feature_data import (DNAIterator, DNAFASTAFormat)


# shamelessly adapted from q2-dada2
_GTE_NEG_1 = (lambda x: x >= -1, 'non-negative; -1 to disable')
_WHOLE_NUM = (lambda x: x >= 0, 'non-negative')
_NAT_NUM = (lambda x: x > 0, 'greater than zero')
_SKIP = (lambda x: True, '')
_valid_inputs = {
    'trim_length': _GTE_NEG_1,
    'mean_error': _NAT_NUM,
    'indel_prob': _NAT_NUM,
    'indel_max': _WHOLE_NUM,
    'min_reads': _WHOLE_NUM,
    'min_size': _WHOLE_NUM,
    'jobs_to_start': _NAT_NUM,
    'hashed_feature_ids': _SKIP,
    'demultiplexed_seqs': _SKIP,
    'reference_seqs': _SKIP
}


# shamelessly adapted from q2-dada2
def _check_inputs(**kwargs):
    for param, arg in kwargs.items():
        check_is_valid, explanation = _valid_inputs[param]
        if not check_is_valid(arg):
            raise ValueError('Argument to %r was %r, should be %s.'
                             % (param, arg, explanation))


def _load_table(base_path):
    """Load the table, remove extraneous filename bits from sample IDs"""
    table = biom.load_table(os.path.join(base_path, 'reference-hit.biom'))
    sid_map = {id_: id_.split('_')[0] for id_ in table.ids(axis='sample')}
    table.update_ids(sid_map, axis='sample', inplace=True)
    return table


def _hash_ids(table):
    """Compute the MD5 of every sequence, update table, return mapping"""
    # Make feature IDs the md5 sums of the sequences.
    # shamelessly adapted from q2-dada2
    fid_map = {id_: hashlib.md5(id_.encode('utf-8')).hexdigest()
               for id_ in table.ids(axis='observation')}
    table.update_ids(fid_map, axis='observation', inplace=True)
    return fid_map


def denoise_16S(
        demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
        trim_length: int,
        mean_error: float=0.005,
        indel_prob: float=0.01,
        indel_max: int=3,
        min_reads: int=10,
        min_size: int=2,
        jobs_to_start: int=1,
        hashed_feature_ids: bool=True) -> (biom.Table, DNAIterator):
    return _denoise_helper(
        demultiplexed_seqs=demultiplexed_seqs,
        mean_error=mean_error,
        indel_prob=indel_prob,
        indel_max=indel_max,
        trim_length=trim_length,
        min_reads=min_reads,
        min_size=min_size,
        jobs_to_start=jobs_to_start,
        hashed_feature_ids=hashed_feature_ids)


def denoise_other(
        demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
        reference_seqs: DNAFASTAFormat,
        trim_length: int,
        mean_error: float=0.005,
        indel_prob: float=0.01,
        indel_max: int=3,
        min_reads: int=10,
        min_size: int=2,
        jobs_to_start: int=1,
        hashed_feature_ids: bool=True) -> (biom.Table, DNAIterator):
    return _denoise_helper(
        demultiplexed_seqs=demultiplexed_seqs,
        reference_seqs=reference_seqs,
        mean_error=mean_error,
        indel_prob=indel_prob,
        indel_max=indel_max,
        trim_length=trim_length,
        min_reads=min_reads,
        min_size=min_size,
        jobs_to_start=jobs_to_start,
        hashed_feature_ids=hashed_feature_ids)


def _denoise_helper(
        demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
        trim_length: int,
        reference_seqs: DNAFASTAFormat=None,
        mean_error: float=0.005,
        indel_prob: float=0.01,
        indel_max: int=3,
        min_reads: int=10,
        min_size: int=2,
        jobs_to_start: int=1,
        hashed_feature_ids: bool=True) -> (biom.Table, DNAIterator):
    _check_inputs(**locals())
    with tempfile.TemporaryDirectory() as tmp:
        seqs_fp = str(demultiplexed_seqs)
        cmd = ['deblur', 'workflow',
               '--seqs-fp', seqs_fp,
               '--output-dir', tmp,
               '--mean-error', str(mean_error),
               '--indel-prob', str(indel_prob),
               '--indel-max', str(indel_max),
               '--trim-length', str(trim_length),
               '--min-reads', str(min_reads),
               '--min-size', str(min_size),
               '-w']

        if reference_seqs is not None:
            cmd.append('--pos-ref-fp')
            cmd.append(str(reference_seqs))

        subprocess.run(cmd, check=True)

        # this is one of the outputs from Deblur, however it isn't clear what
        # the utility of it is for the majority of qiime2 users. on the other
        # hand, it is very easy to test to see if the run completed.
        all_seqs = os.path.join(tmp, 'all.seqs.fa')
        if os.stat(all_seqs).st_size == 0:
            raise ValueError("No sequences passed the filter. It is possible "
                             "the trim_length (%d) may exceed the longest "
                             "sequence, that all of the sequences are "
                             "artifacts like PhiX or adapter, or that the "
                             "positive reference used is not representative "
                             "of the data being denoised.")

        table = _load_table(tmp)

        if hashed_feature_ids:
            obs_map = _hash_ids(table)  # inplace operation
        else:
            obs_map = {i: i for i in table.ids(axis='observation')}

        rep_sequences = DNAIterator(
            (skbio.DNA(k, metadata={'id': v}, lowercase='ignore')
             for k, v in obs_map.items()))

    return (table, rep_sequences)
