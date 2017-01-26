#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    23.01.2017 10:41:07 CET
# File:    test.py

import kdotp_symmetry as kp

if __name__ == '__main__':
    c2 = kp.SymmetryOperation(
        kmatrix=[[0, 1, 0], [1, 0, 0], [0, 0, 1]],
        repr=kp.Representation(
            matrix=[[0, 1], [1, 0]],
            complex_conjugate=False
        )
    )
    print(kp.symmetric_hamiltonian(
        c2,
        expr_basis=kp.monomial_basis(*range(2))
    ))
