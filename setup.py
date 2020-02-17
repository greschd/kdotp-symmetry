# © 2017-2018, ETH Zurich, Institut für Theoretische Physik
# Author:  Dominik Gresch <greschd@gmx.ch>
"""
Usage: pip install -e .[dev]
"""

import re
import sys
if sys.version_info < (3, 6):
    raise 'must use Python version 3.6 or higher'

from setuptools import setup

README = """A tool for calculating the general form of a k.p Hamiltonian under given symmetry constraints."""

with open('./kdotp_symmetry/__init__.py', 'r') as f:
    MATCH_EXPR = "__version__[^'\"]+(['\"])([^'\"]+)"
    VERSION = re.search(MATCH_EXPR, f.read()).group(2)

setup(
    name='kdotp-symmetry',
    version=VERSION,
    url='https://kdotp-symmetry.greschd.ch',
    author='Dominik Gresch',
    author_email='greschd@gmx.ch',
    description=
    'Calculating the general form of a k.p Hamiltonian with given symmetry constraints.',
    install_requires=[
        'sympy<1.5', 'numpy', 'scipy', 'fsc.export',
        'symmetry-representation>=0.3', 'networkx>=2'
    ],
    python_requires=">=3.6",
    extras_require={
        'dev': [
            'pytest', 'yapf==0.29', 'pre-commit', 'prospector==1.2.0',
            'pylint==2.4.4', 'sphinx', 'sphinx_rtd_theme'
        ]
    },
    long_description=README,
    classifiers=[
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English', 'Operating System :: Unix',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'Development Status :: 4 - Beta'
    ],
    license='Apache 2.0',
    packages=['kdotp_symmetry']
)
