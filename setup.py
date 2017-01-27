# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages
import re
import ast

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('q2_deblur/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

setup(
    name="q2-deblur",
    version=version,
    packages=find_packages(),
    install_requires=['qiime >= 2.0.5', 'pandas', 'q2-types >= 0.0.5',
                      'deblur >= 0.1.4'],
    author="Daniel McDonald",
    author_email="mcdonadt@colorado.edu",
    description="Wrapper for Deblur",
    entry_points={
        "qiime.plugins":
        ["q2-deblur=q2_deblur.plugin_setup:plugin"]
    }
)
