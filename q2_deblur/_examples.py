# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

m_p_base = 'https://docs.qiime2.org/{epoch}/data/tutorials/moving-pictures/'
demuxed_seqs_url = m_p_base + 'demux-filtered.qza'

def denoise_16S_example(use):
    demuxed = use.init_artifact_from_url('demux-filtered', demuxed_seqs_url)

    rep_seqs, table, stats = use.action(
        use.UsageAction('deblur', 'denoise_16S'),
        use.UsageInputs(
            demultiplexed_seqs=demuxed,
            trim_length=120,
            sample_stats=True,
        ),
        use.UsageOutputNames(
            representative_sequences='representative_sequences',
            table='table',
            stats='denoising_stats'
        )
    )

    rep_seqs.assert_output_type('FeatureData[Sequence]')
    table.assert_output_type('FeatureTable[Frequency]')
    stats.assert_output_type('SampleData[DeblurStats]')
