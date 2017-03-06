# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model


STATS_HEADER = ['sample-id',

                # number of reads within the per-sample file
                'reads-raw',

                # number of unique reads following dereplication
                # IMPORTANT: this value is impacted by minsize setting. By
                # default, singletons are not included.
                'unique-reads-derep',
                'reads-derep',

                # number of unique reads following deblur
                # IMPORTANT: Deblur corrects frequencies, so the read counts
                # may not fully add up.
                'unique-reads-deblur',
                'reads-deblur',

                # number of reads which appear to be artifact.
                # IMPORTANT: the negative filter is assessed prior to Deblur
                # and by default recruits PhiX and adapter.
                'unique-reads-hit-artifact',
                'reads-hit-artifact',

                # number of reads which appear to be chimeric.
                # IMPORTANT: this is assessed following deblur.
                'unique-reads-chimeric',
                'reads-chimeric',

                # number of reads which appear to be composed of the
                # reference type.
                # IMPORTANT: this is assessed following Deblur, and following
                # the removal of reads below the min-reads threshold.
                'unique-reads-hit-reference',
                'reads-hit-reference',

                # number of reads which do not appear to be composed of
                # the reference type.
                # IMPORTANT: this is assessed following Deblur, and following
                # the removal of reads below the min-reads threshold.
                'unique-reads-missed-reference',
                'reads-missed-reference']


class DeblurStatsFmt(model.TextFileFormat):
    def sniff(self):
        line = open(str(self)).readline()
        hdr = line.strip().split(',')

        return hdr == STATS_HEADER


DeblurStatsDirFmt = model.SingleFileDirectoryFormat(
    'DeblurStatsDirFmt', 'stats.csv', DeblurStatsFmt)
