# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import os
import pkg_resources
import shutil

import pandas as pd
import q2templates

from q2_deblur._format import STATS_DESCRIPTIONS

TEMPLATES = pkg_resources.resource_filename('q2_deblur', 'assets')

COMPUTED_DESCRIPTIONS = {
    'fraction-artifact': ('The fraction of reads which appear to be '
                          'artifactual. This is computed as '
                          'reads-hit-artifact / reads-raw'),
    'fraction-artifact-with-minsize': (
        'The fraction of reads which appear to be artifactual including those '
        'below the min-size threshold. This is computed as '
        '(reads-hit-artifact + (reads-raw - reads-derep)) / reads-raw'),
    'fraction-missed-reference': (
        'The fraction of reads which failed to recruit to the reference. This '
        'is computed as '
        'reads-missed-reference / (reads-deblur - reads-chimeric)')}


def visualize_stats(output_dir: str, deblur_stats: pd.DataFrame) -> None:
    total_artifact = deblur_stats['reads-hit-artifact']
    total_input = deblur_stats['reads-raw']
    deblur_stats['fraction-artifact'] = total_artifact / total_input

    minsize_drop = deblur_stats['reads-raw'] - deblur_stats['reads-derep']
    deblur_stats['fraction-artifact-with-minsize'] = \
        (total_artifact + minsize_drop) / total_input

    total_not_ref = deblur_stats['reads-missed-reference']
    total_deblur = deblur_stats['reads-deblur']
    total_chim = deblur_stats['reads-chimeric']
    deblur_stats['fraction-missed-reference'] = \
        total_not_ref / (total_deblur - total_chim)

    # reorder such that retained fractions follow total-input-reads and
    # total-retained-reads
    columns = list(deblur_stats.columns)[:-3]
    columns.insert(1, 'fraction-missed-reference')
    columns.insert(1, 'fraction-artifact')
    columns.insert(1, 'fraction-artifact-with-minsize')
    deblur_stats = deblur_stats[columns]

    deblur_stats.sort_values('fraction-artifact-with-minsize', inplace=True,
                             ascending=False)

    deblur_stats = deblur_stats.reset_index()
    html = deblur_stats.to_html(classes='table table-striped table-hover')
    html = html.replace('border="1"', 'border="0"')
    html = html.replace('table-hover"', 'table-hover" id="stats"')

    # ghetto force in tooltips
    description_sources = STATS_DESCRIPTIONS.copy()
    description_sources.update(COMPUTED_DESCRIPTIONS)
    htmlparts = html.splitlines()
    headstart = None
    headend = None
    for idx, line in enumerate(htmlparts):
        if '<thead>' in line:
            headstart = idx
        elif '</thead>' in line:
            headend = idx

    regex = re.compile("<th>(.*?)</th>")
    new_header = []
    for entry in htmlparts[headstart:headend]:
        new_entry = entry[:]
        if '<th>' in entry and entry.strip() != '<th></th>':
            label = regex.findall(entry)[0]
            desc = description_sources[label]
            label = '<th data-toggle="tooltip" title="%s">%s</th>' % (desc,
                                                                      label)
            new_entry = label
        new_header.append(new_entry)
    htmlparts[headstart:headend] = new_header
    html = '\n'.join(htmlparts)

    index = os.path.join(TEMPLATES, 'index.html')

    context = {
        'result': html
    }

    js = os.path.join(TEMPLATES, 'js', 'tsorter.min.js')
    os.mkdir(os.path.join(output_dir, 'js'))
    shutil.copy(js, os.path.join(output_dir, 'js', 'tsorter.min.js'))

    q2templates.render(index, output_dir, context=context)
