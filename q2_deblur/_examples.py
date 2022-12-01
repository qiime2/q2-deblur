# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


demuxed_seqs_url = ('https://data.qiime2.org/usage-examples/'
                    'moving-pictures/demux-filtered.qza')
denoise_stats_url = ('https://data.qiime2.org/usage-examples/'
                     'moving-pictures/deblur-stats.qza')


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
    stats.assert_output_type('DeblurStats')


def visualize_stats_example(use):
    stats = use.init_artifact_from_url('deblur-stats', denoise_stats_url)

    viz, = use.action(
        use.UsageAction('deblur', 'visualize_stats'),
        use.UsageInputs(
            deblur_stats=stats
        ),
        use.UsageOutputNames(
            visualization='deblur-stats'
        )
    )

    viz.assert_output_type('Visualization')
