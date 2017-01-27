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
    version="2017.1.0.dev0",
    packages=find_packages(),
    install_requires=['qiime >= 2.0.6', 'pandas', 'q2-types >= 2.0.6',
                      'deblur >= 0.1.8'],
    author="Daniel McDonald",
    author_email="mcdonadt@colorado.edu",
    description="Wrapper for Deblur",
    entry_points={
        "qiime.plugins":
        ["q2-deblur=q2_deblur.plugin_setup:plugin"]
    }
)
