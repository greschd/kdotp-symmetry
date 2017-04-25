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

with open('./kdotp_symmetry/_version.py', 'r') as f:
    match_expr = "__version__[^'" + '"]+([' + "'" + r'"])([^\1]+)\1'
    version = re.search(match_expr, f.read()).group(2)

setup(
    name='kdotp-symmetry',
    version=version,
    url='http://z2pack.ethz.ch/kdotp-symmetry',
    author='Dominik Gresch',
    author_email='greschd@gmx.ch',
    description='Calculating the general form of a k.p Hamiltonian with given symmetry constraints.',
    install_requires=['sympy', 'numpy', 'scipy', 'fsc.export', 'symmetry-representation'],
    long_description=readme,
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'Development Status :: 4 - Beta'
    ],
    license='GPL',
    packages=['kdotp_symmetry']
)
