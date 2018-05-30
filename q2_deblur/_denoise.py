# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import tempfile
import subprocess
import os
import hashlib
import gzip
import itertools

import numpy as np
import biom
import skbio
import pandas as pd
from q2_types.per_sample_sequences import \
        SingleLanePerSampleSingleEndFastqDirFmt, FastqGzFormat
from q2_types.feature_data import (DNAIterator, DNAFASTAFormat)

from q2_deblur._format import STATS_HEADER


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
    'reference_seqs': _SKIP,
    'sample_stats': _SKIP
}


# shamelessly adapted from q2-dada2
def _check_inputs(**kwargs):
    for param, arg in kwargs.items():
        check_is_valid, explanation = _valid_inputs[param]
        if not check_is_valid(arg):
            raise ValueError('Argument to %r was %r, should be %s.'
                             % (param, arg, explanation))


def _load_table(base_path, name='reference-hit.biom'):
    """Load the table, remove extraneous filename bits from sample IDs"""
    table = biom.load_table(os.path.join(base_path, name))
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
        sample_stats: bool=False,
        mean_error: float=0.005,
        indel_prob: float=0.01,
        indel_max: int=3,
        min_reads: int=10,
        min_size: int=2,
        jobs_to_start: int=1,
        hashed_feature_ids: bool=True) -> (biom.Table,
                                           DNAIterator,
                                           pd.DataFrame):
    return _denoise_helper(
        sample_stats=sample_stats,
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
        sample_stats: bool=False,
        mean_error: float=0.005,
        indel_prob: float=0.01,
        indel_max: int=3,
        min_reads: int=10,
        min_size: int=2,
        jobs_to_start: int=1,
        hashed_feature_ids: bool=True) -> (biom.Table,
                                           DNAIterator,
                                           pd.DataFrame):
    return _denoise_helper(
        sample_stats=sample_stats,
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
        sample_stats: bool=False,
        reference_seqs: DNAFASTAFormat=None,
        mean_error: float=0.005,
        indel_prob: float=0.01,
        indel_max: int=3,
        min_reads: int=10,
        min_size: int=2,
        jobs_to_start: int=1,
        hashed_feature_ids: bool=True) -> (biom.Table,
                                           DNAIterator,
                                           pd.DataFrame):
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
               '--jobs-to-start', str(jobs_to_start),
               '-w']

        if reference_seqs is not None:
            cmd.append('--pos-ref-fp')
            cmd.append(str(reference_seqs))

        if sample_stats:
            cmd.append('--keep-tmp-files')

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

        if sample_stats:
            stats = _gather_stats(demultiplexed_seqs, tmp)
        else:
            stats = pd.DataFrame([], columns=STATS_HEADER)
            stats.set_index('sample-id', inplace=True)

    return (table, rep_sequences, stats)


def _gather_stats(demux, tmp):
    workingdir = os.path.join(tmp, 'deblur_working_dir')
    if not os.path.exists(workingdir):
        raise IOError("Cannot find the deblur_working_dir")

    demux_manifest = demux.manifest.view(demux.manifest.format)
    demux_manifest = pd.read_csv(demux_manifest.open(), dtype=str)
    demux_manifest.set_index('filename', inplace=True)

    all_table = _load_table(tmp, 'all.biom')
    ref_table = _load_table(tmp, 'reference-hit.biom')

    stats = []

    iterator = demux.sequences.iter_views(FastqGzFormat)
    for bc_id, (fname, fp) in enumerate(iterator):
        sample_id = demux_manifest.loc[str(fname)]['sample-id']

        # actual raw read counts
        raw_counts = _read_fastq_seqs(str(fp))

        # VSEARCH dereplicated raw input, sum of the reads represented
        unique_reads_derep, reads_derep = _fasta_counts(workingdir, fname,
                                                        'trim.derep')

        # VSEARCH dereplicated raw input minus sequences which
        # recruit to the negative reference database. By default
        # the negative reference database is composed of PhiX and
        # a 16S Illumina adapter.
        nonartifact_unique, nonartifact_counts = \
            _fasta_counts(workingdir, fname, 'trim.derep.no_artifacts')

        reads_hit_artifact = reads_derep - nonartifact_counts
        unique_reads_hit_artifact = unique_reads_derep - nonartifact_unique

        unique_reads_deblur, reads_deblur = \
            _fasta_counts(workingdir, fname,
                          'trim.derep.no_artifacts.msa.deblur')

        unique_reads_chim, reads_chim = \
            _fasta_counts(workingdir, fname,
                          'trim.derep.no_artifacts.msa.deblur.no_chimeras')
        unique_reads_chim = unique_reads_deblur - unique_reads_chim
        reads_chim = reads_deblur - reads_chim

        if all_table.exists(sample_id):
            all_data = all_table.data(sample_id, dense=False)
        else:
            all_data = np.zeros((1))

        if ref_table.exists(sample_id):
            ref_data = ref_table.data(sample_id, dense=False)
        else:
            ref_data = np.zeros((1))

        final_reads_deblur = all_data.sum()
        final_unique_reads_deblur = (all_data > 0).sum()

        reads_hit_ref = ref_data.sum()
        unique_reads_hit_ref = (ref_data > 0).sum()

        reads_missed_ref = final_reads_deblur - reads_hit_ref
        unique_reads_missed_ref = \
            final_unique_reads_deblur - unique_reads_hit_ref

        stats.append((sample_id, raw_counts,
                      unique_reads_derep, reads_derep,
                      unique_reads_deblur, reads_deblur,
                      unique_reads_hit_artifact, reads_hit_artifact,
                      unique_reads_chim, reads_chim,
                      unique_reads_hit_ref, reads_hit_ref,
                      unique_reads_missed_ref, reads_missed_ref))

    df = pd.DataFrame(stats, columns=STATS_HEADER, dtype=object)
    df = df.set_index('sample-id')
    return df.astype(int)


def _fasta_counts(workingdir, sample_id, suffix):
    # This function is adapted from @jairideout's SO post:
    # http://stackoverflow.com/a/39302117/3424666
    counts = 0
    unique = 0

    path = os.path.join(workingdir, '%s.%s' % (sample_id, suffix))
    if not os.path.exists(path):
        return 0, 0

    with open(os.path.join(path)) as fh:
        for seq_header, seq in itertools.zip_longest(*[fh] * 2):
            # >foo stuff;size=123;
            size = seq_header.rsplit(';', 2)[1]
            counts += int(size.split('=')[1])
            unique += 1
        return unique, counts


def _read_fastq_seqs(filepath):
    if not os.path.exists(filepath):
        return 0

    fh = gzip.open(filepath, 'rt')
    return sum(1 for _ in fh) / 4
