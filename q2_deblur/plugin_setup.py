# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess
import tempfile
import os
import qiime2.plugin
import biom
import skbio
import hashlib

from q2_types.feature_table import FeatureTable, Frequency
from q2_types.per_sample_sequences import \
        SingleLanePerSampleSingleEndFastqDirFmt
from q2_types.per_sample_sequences import SequencesWithQuality
from q2_types.feature_data import DNAIterator
from q2_types.sample_data import SampleData
from q2_types.feature_data import FeatureData, Sequence
import q2_deblur


def denoise(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
            pos_ref_filepath: str=None,
            neg_ref_filepath: str=None,
            mean_error: float=0.005,
            indel_prob: float=0.01,
            indel_max: int=3,
            trim_length: int=150,
            min_reads: int=0,
            min_size: int=2,
            negate: bool=False,
            jobs_to_start: int=1,
            hashed_feature_ids: bool=True) -> (biom.Table, DNAIterator):

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
        if pos_ref_filepath is not None:
            cmd.append('--pos-ref-db')
            cmd.append(pos_ref_filepath)

        if neg_ref_filepath is not None:
            cmd.append('--neg-ref-db')
            cmd.append(neg_ref_filepath)

        if negate:
            cmd.append('--negate')

        subprocess.run(cmd, check=True)

        # code adapted from q2-dada2
        table = biom.load_table(os.path.join(tmp, 'final.biom'))
        sid_map = {id_: id_.split('_')[0] for id_ in table.ids(axis='sample')}
        table.update_ids(sid_map, axis='sample', inplace=True)

        if hashed_feature_ids:
            # Make feature IDs the md5 sums of the sequences.
            fid_map = {id_: hashlib.md5(id_.encode('utf-8')).hexdigest()
                       for id_ in table.ids(axis='observation')}
            table.update_ids(fid_map, axis='observation', inplace=True)

            rep_sequences = DNAIterator((skbio.DNA(k, metadata={'id': v},
                                                   lowercase='ignore')
                                         for k, v in fid_map.items()))
        else:
            rep_sequences = DNAIterator(
                (skbio.DNA(id_, metadata={'id': id_}, lowercase='ignore')
                 for id_ in table.ids(axis='observation')))

    return (table, rep_sequences)


plugin = qiime2.plugin.Plugin(
    name='deblur',
    version=q2_deblur.__version__,
    website='https://github.com/biocore/deblur',
    package='q2_deblur',
    user_support_text='https://github.com/biocore/deblur/issues',
    citation_text=("Deblur rapidly resolves single-nucleotide community "
                   "sequence patterns. Amnon Amir, Daniel McDonald, Jose "
                   "A. Navas-Molina, Evguenia Kopylova, Jamie Morton, "
                   "Zhenjiang Zech Xu, Eric P. Kightley, Luke R. Thompson, "
                   "Embriette R. Hyde, Antonio Gonzalez, Rob Knight. mSystems "
                   "(in press).")
)

plugin.methods.register_function(
    function=denoise,
    inputs={
        'demultiplexed_seqs': SampleData[SequencesWithQuality]
    },
    parameters={
        'pos_ref_filepath': qiime2.plugin.Str,
        'neg_ref_filepath': qiime2.plugin.Str,
        'mean_error': qiime2.plugin.Float,
        'indel_prob': qiime2.plugin.Float,
        'indel_max': qiime2.plugin.Int,
        'trim_length': qiime2.plugin.Int,
        'min_reads': qiime2.plugin.Int,
        'min_size': qiime2.plugin.Int,
        'negate': qiime2.plugin.Bool,
        'jobs_to_start': qiime2.plugin.Int,
        'hashed_feature_ids': qiime2.plugin.Bool
    },
    parameter_descriptions={
        'pos_ref_filepath': ("Positive filtering database. Keep all sequences "
                             "aligning to these sequences."),
        'neg_ref_filepath': ("Negative (artifacts) filtering database. "
                             "Discard sequences aligning to these sequences."),
        'mean_error': ("The mean per nucleotide error, used for original "
                       "sequence estimate. If not passed, a value of 0.5% is "
                       "used."),
        'indel_prob': ('Insertion/deletion (indel) probability (same for N '
                       'indels).'),
        'indel_max': "Maximum number of insertion/deletions.",
        'trim_length': "Sequence trim length.",
        'min_reads': ("Retain only features appearing at least min_reads "
                      "times across all samples in the resulting feature "
                      "table."),
        'min_size': ("In each sample, discard all features with an abundance "
                     "less than min_size."),
        'negate': ("Discard all sequences aligning to "
                   "the sequences provided in neg_ref_fp. "
                   "Used for removal of phiX and adapter sequences."
                   " Note an additional positive filtering file "
                   "(only sequences close enough to pos_ref_filepath) will "
                   "also be generated from these reads."),
        'jobs_to_start': "Number of jobs to start (if to run in parallel).",
        'hashed_feature_ids': "If true, hash the feature IDs."
    },
    outputs=[('table', FeatureTable[Frequency]),
             ('representative_sequences', FeatureData[Sequence])],
    output_descriptions={
        'table': 'The resulting denoised feature table.',
        'representative_sequences': 'The resulting feature sequences.'
    },
    name='Perform sequence quality control using the deblur workflow',
    description='Apply the deblur quality control workflow'
)
