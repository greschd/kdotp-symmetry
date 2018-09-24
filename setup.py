"""
pip install -e .[dev]
"""

import re
import sys
if sys.version_info < (3, 5):
    raise 'must use Python version 3.5 or higher'

from setuptools import setup

README = """A tool for calculating the general form of a k.p Hamiltonian under given symmetry constraints."""

with open('./kdotp_symmetry/__init__.py', 'r') as f:
    MATCH_EXPR = "__version__[^'\"]+(['\"])([^'\"]+)"
    VERSION = re.search(MATCH_EXPR, f.read()).group(2)

setup(
    name='kdotp-symmetry',
    version=VERSION,
    url='http://z2pack.ethz.ch/kdotp-symmetry',
    author='Dominik Gresch',
    author_email='greschd@gmx.ch',
    description=
    'Calculating the general form of a k.p Hamiltonian with given symmetry constraints.',
    install_requires=[
        'sympy', 'numpy', 'scipy', 'fsc.export', 'symmetry-representation>=0.1'
    ],
    extras_require={'dev': ['pytest', 'yapf', 'pre-commit', 'prospector']},
    long_description=README,
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English', 'Operating System :: Unix',
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
