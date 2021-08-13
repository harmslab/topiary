#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import sys
from shutil import rmtree

from setuptools import find_packages, setup, Command, Extension

# Package meta-data.
NAME = 'topiary'
DESCRIPTION = 'A lightweight python package for doing phylogenetics work (tailored for ancestral sequence reconstruction).'
URL = 'https://github.com/harmslab/topiary'
EMAIL = 'harmsm@gmail.com'
AUTHOR = 'Mike Harms'
REQUIRES_PYTHON = '>=3.6.0'
VERSION = None

# What packages are required for this module to be executed?
REQUIRED = [
    "biopython>=1.79",
    "ete3>=3.1.2",
    "opentree>=1.0.1",
    "tqdm>=4.61.2",
    "dendropy>=4.5.2",
    "numpy>=1.21.1",
    "pandas>=1.3.1",
    "matplotlib>=3.4.2",
    "pastml>=1.9.34",
]


here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = '\n' + f.read()

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
    with open(os.path.join(here, NAME, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION

# Where the magic happens:
setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=('tests',)),
    install_requires=REQUIRED,
    extras_require = {
        'test': ['pytest'],
    },
    scripts=['bin/run-raxml'],
    include_package_data=True,
    license='MIT',
    classifiers=[
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Science/Research',
      'Programming Language :: Python :: 3.6',
      'Programming Language :: Python :: 3.7',
      'Programming Language :: Python :: 3.8',
      'Programming Language :: Python :: 3.9',
    ],
    keywords='phylogenetics ASR'
)
