#!/usr/bin/env python

import sys
from setuptools import setup

"""setup.py: setuptools control."""

install_requires = [
    'numpy>=1.7.0',
    'netCDF4>=1.0.8',
]

if sys.version_info < (2, 7):
    install_requires.append('argparse>=1.2.1')

setup(
    name='layerpack',
    version='0.1.0',
    description='Command-line tools for packing/unpacking NetCDF arrays in slices',
    long_description=open('README.rst').read(),
    author='Jeremy Silver',
    author_email='jeremy.silver@unimelb.edu.au',
##    url='http://layerpack.rtfd.org/',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Utilities',
        'Topic :: System :: Archiving :: Compression'
    ],
    packages=[
        'layerpack'
    ],
    entry_points={
        'console_scripts': [
            'ncpack = layerpack.ncpack:main',
            'ncunpack = layerpack.ncunpack:main',
            'nccheckdiff = layerpack.nccheckdiff:main'
        ]
    },
    install_requires=install_requires
)
