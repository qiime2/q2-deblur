# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from setuptools import setup, find_packages

setup(
    name="q2-deblur",
    version="2017.2.3",
    packages=find_packages(),
    install_requires=['qiime2 == 2017.2.*', 'pandas', 'q2-types == 2017.2.*',
                      'deblur >= 1.0.1'],
    author="Daniel McDonald",
    author_email="wasade@gmail.com",
    description="Sequence quality control with deblur",
    entry_points={
        "qiime2.plugins":
        ["q2-deblur=q2_deblur.plugin_setup:plugin"]
    }
)
