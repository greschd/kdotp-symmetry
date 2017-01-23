#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    20.10.2014 11:27:40 CEST
# File:    setup.py

import re
try:
    from setuptools import setup
except:
    from distutils.core import setup

import sys
if sys.version_info < (3, 5):
    raise 'must use Python version 3.5 or higher'

readme = """A tool for calculating the general form of a k.p Hamiltonian under given symmetry constraints."""

with open('./tbmodels/_version.py', 'r') as f:
    match_expr = "__version__[^'" + '"]+([' + "'" + r'"])([^\1]+)\1'
    version = re.search(match_expr, f.read()).group(2)

setup(
    name='tbmodels',
    version=version,
    author='Dominik Gresch',
    author_email='greschd@gmx.ch',
    description='Calculating the general form of a k.p Hamiltonian with given symmetry constraints.',
    install_requires=['sympy'],
    long_description=readme,
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics'
    ],
    license='GPL',
    packages=['kdotp_symmetry']
)
