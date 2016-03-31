#!/usr/bin/env python
import relic.release
from glob import glob
from numpy import get_include as np_include
from setuptools import setup, find_packages, Extension


version = relic.release.get_info()
relic.release.write_template(version, 'lib/wfpc2tools')

setup(
    name = 'wfpc2tools',
    version = version.pep386,
    author = 'Warren Hack, Nadezhda Dencheva, David Grumm',
    author_email = 'help@stsci.edu',
    description = 'Tools for use with WFPC2 (Wide Field and Planetary Camera 2)',
    url = 'https://github.com/spacetelescope/wfpc2tools',
    classifiers = [
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    install_requires = [
        'astropy',
        'nose',
        'numpy',
        'scipy',
        'sphinx',
        'stsci.imagestats',
        'stsci.sphinxext',
        'stsci.tools',
    ],
    package_dir = {
        '': 'lib',
    },
    packages = find_packages('lib'),
)
