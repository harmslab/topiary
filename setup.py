#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import glob

from setuptools import find_packages, setup

# Package meta-data.
DESCRIPTION = \
"""Python framework for doing ancestral sequence reconstruction."""
URL = "https://github.com/harmslab/topiary"
EMAIL = "harmsm@gmail.com"
AUTHOR = "Michael J. Harms"
REQUIRES_PYTHON = ">=3.8.0"
VERSION = None

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = '\n' + f.read()

# Load the package's __version__.py module as a dictionary.
about = {}
if VERSION is None:
    with open(os.path.join(here, "topiary", '__version__.py')) as f:
        exec(f.read(),about)
else:
    about['__version__'] = VERSION

# Where the magic happens:
setup(
    name="topiary-asr", # Note this is different than "topiary" for PIP package name
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=('tests',)),
    scripts=glob.glob("bin/topiary-*"),
    include_package_data=True,
    license='MIT',
    classifiers=[
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Science/Research',
      'Programming Language :: Python :: 3.8',
      'Programming Language :: Python :: 3.9',
      'Programming Language :: Python :: 3.10',
    ],
    keywords="phylogenetics; ancestral sequence reconstruction; ASR; bioinformatics; protein; evolutionary biochemistry"
)
