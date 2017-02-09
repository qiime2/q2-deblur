# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import qiime2.plugin
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Sequence
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (SequencesWithQuality,
                                           PairedEndSequencesWithQuality)

import q2_deblur


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


_parameter_descriptions = {
    'mean_error': ("The mean per nucleotide error, used for original "
                   "sequence estimate."),
    'indel_prob': ('Insertion/deletion (indel) probability (same for N '
                   'indels).'),
    'indel_max': "Maximum number of insertion/deletions.",
    'trim_length': "Sequence trim length, specify -1 to disable trimming.",
    'min_reads': ("Retain only features appearing at least min_reads "
                  "times across all samples in the resulting feature "
                  "table."),
    'min_size': ("In each sample, discard all features with an abundance "
                 "less than min_size."),
    'jobs_to_start': "Number of jobs to start (if to run in parallel).",
    'hashed_feature_ids': "If true, hash the feature IDs."
}


_output_descriptions = {
    'table': 'The resulting denoised feature table.',
    'representative_sequences': 'The resulting feature sequences.'
}


_parameters = {
    'mean_error': qiime2.plugin.Float,
    'indel_prob': qiime2.plugin.Float,
    'indel_max': qiime2.plugin.Int,
    'trim_length': qiime2.plugin.Int,
    'min_reads': qiime2.plugin.Int,
    'min_size': qiime2.plugin.Int,
    'jobs_to_start': qiime2.plugin.Int,
    'hashed_feature_ids': qiime2.plugin.Bool
}


_outputs = [('table', FeatureTable[Frequency]),
            ('representative_sequences', FeatureData[Sequence])]


plugin.methods.register_function(
    function=q2_deblur.denoise_16S,
    inputs={
        'demultiplexed_seqs': SampleData[SequencesWithQuality |
                                         PairedEndSequencesWithQuality],
    },
    parameters=_parameters,
    outputs=_outputs,
    input_descriptions={
        'demultiplexed_seqs': 'The demultiplexed sequences to be denoised.',
    },
    parameter_descriptions=_parameter_descriptions,
    output_descriptions=_output_descriptions,
    name='Deblur sequences using a 16S positive filter.',
    description=('Perform sequence quality control for Illumina data using '
                 'the Deblur workflow with a 16S reference as a positive '
                 'filter. Only forward reads are supported at this time. The '
                 'specific reference used is the 88% OTUs from Greengenes '
                 '13_8. This mode of operation should only be used when data '
                 'were generated from a 16S amplicon protocol on an Illumina '
                 'platform. The reference is only used to assess whether each '
                 'sequence is likely to be 16S by a local alignment using '
                 'SortMeRNA with a permissive e-value; the reference is not '
                 'used to characterize the sequences.')
)


plugin.methods.register_function(
    function=q2_deblur.denoise_other,
    inputs={
        'demultiplexed_seqs': SampleData[SequencesWithQuality |
                                         PairedEndSequencesWithQuality],
        'reference_seqs': FeatureData[Sequence],
    },
    parameters=_parameters,
    outputs=_outputs,
    input_descriptions={
        'demultiplexed_seqs': 'The demultiplexed sequences to be denoised.',
        'reference_seqs': ("Positive filtering database. Keep all "
                           "sequences aligning to these sequences."),
    },
    parameter_descriptions=_parameter_descriptions,
    output_descriptions=_output_descriptions,
    name='Deblur sequences using a user-specified positive filter.',
    description=('Perform sequence quality control for Illumina data using '
                 'the Deblur workflow, including positive alignment-based '
                 'filtering. Only forward reads are supported at this time. '
                 'This mode of execution is particularly useful when '
                 'operating on non-16S data. For example, to apply Deblur to '
                 '18S data, you would want to specify a reference composed of '
                 '18S sequences in order to filter out sequences which do not '
                 'appear to be 18S. The assessment is performed by local '
                 'alignment using SortMeRNA with a permissive e-value '
                 'threshold.')
)
