#!/usr/bin/env python

from setuptools import find_packages
from setuptools import setup

PACKAGENAME = 'wfpc2tools'

setup(
    name=PACKAGENAME,
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    install_requires=[
        'astropy>=0.3.1',
        'numpy>=1.5.1',
        'scipy>=0.14',
        'stsci.imagestats',
        'stsci.tools',
    ],
    extras_require={
        'docs': [
            'sphinx',
        ],
        'test': [
            'pytest',
            'pytest-cov',
            'codecov',
        ],
    },
    packages=find_packages(),
    author='Warren Hack, Nadezhda Dencheva, David Grumm',
    author_email='help@stsci.edu',
    description='Tools for use with WFPC2 (Wide Field and Planetary Camera 2)',
    url='https://github.com/spacetelescope/wfpc2tools',
    classifiers = [
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)
