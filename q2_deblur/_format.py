# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
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

STATS_DESCRIPTIONS = {
    'sample-id': "The sample ID",
    'reads-raw': "The number of reads presented to Deblur",
    'unique-reads-derep': "The number of unique reads following dereplicaton",
    'reads-derep': ("The number of reads following dereplication. IMPORTANT: "
                    "the minsize parameter is accounted at this step; if using"
                    " defaults, singletons will not be included here."),
    'unique-reads-deblur': "The number of unique reads following Deblur",
    'reads-deblur': ("The number of reads following Deblur. IMPORTANT: Deblur "
                     "adjusts frequencies so the read counts after factoring "
                     "in singletons may not sum to the total input reads"),
    'unique-reads-hit-artifact': ("The number of unique reads which recruited "
                                  "to the negative filter database. This is "
                                  "assessed following the application of "
                                  "Deblur"),
    'reads-hit-artifact': ("The number of reads which recruited to the "
                           "negative filter database. This is assessed "
                           "following the application of Deblur."),
    'unique-reads-chimeric': ("The number of unique reads which appear to be "
                              "chimeric. This is assessed following Deblur."),
    'reads-chimeric': ("The number of reads which appear to be chimeric. This "
                       "is assessed following Deblur."),
    'unique-reads-hit-reference': ("The number of unique Deblur reads which "
                                   "recruited to the positive reference. "
                                   "IMPORTANT: this is assessed after the "
                                   "removal of features below the min-reads "
                                   "threshold."),
    'reads-hit-reference': ("The number of Deblur reads which recruited to the"
                            " positive reference. IMPORTANT: this is "
                            "assessed after the removal of features below the "
                            "min-reads threshold."),
    'unique-reads-missed-reference': ("The number of unique Deblur reads "
                                      "which failed to recruit to the positive"
                                      " reference. IMPORTANT: this is "
                                      "assessed after the removal of features "
                                      "below the min-reads threshold."),
    'reads-missed-reference': ("The number of Deblur reads which failed to "
                               "recruit to the positive reference. "
                               "IMPORTANT: this is assessed after the "
                               "removal of features below the min-reads "
                               "threshold.")
    }


class DeblurStatsFmt(model.TextFileFormat):
    def sniff(self):
        line = open(str(self)).readline()
        hdr = line.strip().split(',')

        return hdr == STATS_HEADER


DeblurStatsDirFmt = model.SingleFileDirectoryFormat(
    'DeblurStatsDirFmt', 'stats.csv', DeblurStatsFmt)
