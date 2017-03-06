# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources

import pandas as pd
import q2templates

TEMPLATES = pkg_resources.resource_filename('q2_deblur', 'assets')


def visualize_stats(output_dir: str, filter_stats: pd.DataFrame) -> None:
    sums = filter_stats.sum()
    sums.name = 'Totals'
    filter_stats = filter_stats.append(sums)

    filter_stats.sort_values('input-reads-raw', inplace=True,
                             ascending=False)

    total_artifact = filter_stats['reads-hit-artifact']
    total_input = filter_stats['input-reads-raw']
    filter_stats['fraction-artifact'] = total_artifact / total_input

    total_artifact = filter_stats['reads-hit-artifact']
    total_artifact += total_input - filter_stats['reads-derep']
    filter_stats['fraction-artifact-with-minsize'] = \
            total_artifact / total_input

    total_ref = filter_stats['reads-hit-reference']
    total_not_ref = filter_stats['reads-missed-reference']
    total_deblur = filter_stats['reads-deblur']
    filter_stats['fraction-missed-reference'] = total_not_ref / total_deblur

    # reorder such that retained fractions follow total-input-reads and
    # total-retained-reads
    columns = list(filter_stats.columns)[:-3]
    columns.insert(1, 'fraction-missed-reference')
    columns.insert(1, 'fraction-artifact-with-minsize')
    columns.insert(1, 'fraction-artifact')
    filter_stats = filter_stats[columns]

    html = filter_stats.to_html(classes='table table-striped table-hover')
    html = html.replace('border="1"', 'border="0"')
    index = os.path.join(TEMPLATES, 'index.html')
    context = {
        'result': html
    }

    q2templates.render(index, output_dir, context=context)

