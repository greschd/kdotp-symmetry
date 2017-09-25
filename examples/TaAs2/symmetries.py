#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    26.01.2017 22:40:33 CET
# File:    symmetries.py

import sympy as sp
from sympy.core.numbers import I
import sympy.physics.matrices as sm
from sympy.physics.quantum import TensorProduct
import symmetry_representation as sr

import kdotp_symmetry as kp

# In this project we used the basis of tensor products of Pauli matrices
pauli_vec = [sp.eye(2), *(sm.msigma(i) for i in range(1, 4))]
basis = [TensorProduct(p1, p2) for p1 in pauli_vec for p2 in pauli_vec]

# creating the symmetry operations
c2y = sr.SymmetryOperation(
    rotation_matrix=[[0, 1, 0], [1, 0, 0], [0, 0, -1]],
    repr_matrix=sp.diag(I, -I, I, -I),
    repr_has_cc=False
)

parity = sr.SymmetryOperation(
    rotation_matrix=-sp.eye(3),
    repr_matrix=sp.diag(1, 1, -1, -1),
    repr_has_cc=False
)

time_reversal = sr.SymmetryOperation(
    rotation_matrix=sp.eye(3),
    repr_matrix=TensorProduct(sp.eye(2), sp.Matrix([[0, -1], [1, 0]])),
    repr_has_cc=True
)


def print_result(order):
    """prints the basis for a given order of k"""
    print('Order:', order)
    for m in kp.symmetric_hamiltonian(
        c2y,
        parity,
        time_reversal,
        expr_basis=kp.monomial_basis(order),
        repr_basis=basis
    ):
        print(m)
    print()


if __name__ == '__main__':
    for i in range(3):
        print_result(order=i)
