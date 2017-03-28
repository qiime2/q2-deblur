# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

import versioneer


setup(
    name="q2-deblur",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    packages=find_packages(),
    author="Daniel McDonald",
    author_email="wasade@gmail.com",
    description="Sequence quality control with deblur",
    entry_points={
        "qiime2.plugins":
        ["q2-deblur=q2_deblur.plugin_setup:plugin"]
    },
    package_data={
        "q2_deblur": ["assets/index.html"]
    },
    zip_safe=True,
)
