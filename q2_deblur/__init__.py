# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._denoise import denoise_16S, denoise_other
from ._viz_stats import visualize_stats

from ._version import get_versions


__version__ = get_versions()['version']
del get_versions

__all__ = ['denoise_16S', 'denoise_other', 'visualize_stats']
