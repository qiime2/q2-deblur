# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pkg_resources

from ._denoise import denoise_16S, denoise_other

__version__ = pkg_resources.get_distribution('q2-deblur').version
__all__ = ['denoise_16S', 'denoise_other']
