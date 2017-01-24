from setuptools import setup, find_packages

setup(
    name="q2-deblur",
    version="2017.2.0.dev0",
    packages=find_packages(),
    install_requires=['qiime2 == 2017.2.*', 'pandas', 'q2-types == 2017.2.*',
                      'deblur >= 0.1.8'],
    author="Daniel McDonald",
    author_email="mcdonadt@colorado.edu",
    description="Sequence quality control with deblur",
    entry_points={
        "qiime2.plugins":
        ["q2-deblur=q2_deblur.plugin_setup:plugin"]
    }
)
