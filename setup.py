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
    # pandas and q2-dummy-types are only required for the dummy methods and
    # visualizers provided as examples. Remove these dependencies when you're
    # ready to develop your plugin, and add your own dependencies (if there are
    # any).
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
