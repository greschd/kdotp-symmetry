#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    23.01.2017 10:41:07 CET
# File:    test.py

from kdotp_symmetry.analyse import *

if __name__ == '__main__':
    c2 = SymmetryOperation(
        kmatrix=[[0, 1, 0], [1, 0, 0], [0, 0, 1]],
        repr=Representation(
            matrix=[[0, 1], [1, 0]],
            complex_conjugate=False
        )
    )
    print(symmetric_hamiltonian(c2, power=2))
