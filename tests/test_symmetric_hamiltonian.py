#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    26.01.2017 14:34:15 CET
# File:    test_symmetric_hamiltonian.py

import pytest
import sympy as sp
from sympy import Matrix
from sympy.core.numbers import I

import kdotp_symmetry as kp

kx, ky, kz = sp.symbols('kx, ky, kz')

@pytest.mark.parametrize('symmetry_operations,expr_basis,repr_basis,result', [
    (
        [kp.SymmetryOperation(
            kmatrix=[[0, 1, 0], [1, 0, 0], [0, 0, 1]],
            repr=kp.Representation(
                matrix=[[0, 1], [1, 0]],
                complex_conjugate=False
            )
        )],
        kp.monomial_basis(0),
        'auto',
        (
            [[1, 1, 0, 0], [0, 0, 1, 0]], 
            [
                Matrix([[1, 0], [0, 0]]), 
                Matrix([[0, 0], [0, 1]]), 
                Matrix([[0, 1], [1, 0]]), 
                Matrix([[0, -I], [I,  0]])
            ], 
            [
                Matrix([[1, 0], [0, 1]]), 
                Matrix([[0, 1], [1, 0]])
            ]
        )
    ),
    (
        [kp.SymmetryOperation(
            kmatrix=[[0, 1, 0], [1, 0, 0], [0, 0, 1]],
            repr=kp.Representation(
                matrix=[[0, 1], [1, 0]],
                complex_conjugate=False
            )
        )],
        kp.monomial_basis(1),
        'auto',
        (
            [[0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]], 
            [
                Matrix([[kx, 0],[ 0, 0]]),
                Matrix([[0,  0],[0, kx]]),
                Matrix([[0, kx],[kx,  0]]),
                Matrix([[0, -I*kx],[I*kx, 0]]),
                Matrix([[ky, 0],[0, 0]]),
                Matrix([[0,  0],[0, ky]]),
                Matrix([[0, ky],[ky, 0]]),
                Matrix([[0, -I*ky],[I*ky, 0]]),
                Matrix([[kz, 0],[0, 0]]),
                Matrix([[0,  0],[0, kz]]),
                Matrix([[0, kz],[kz, 0]]),
                Matrix([[0, -I*kz],[I*kz, 0]])
            ], 
            [
                Matrix([[ky, 0],[0, kx]]),
                Matrix([[kx,  0],[0, ky]]),
                Matrix([[0, kx + ky],[kx + ky, 0]]),
                Matrix([[0, I*kx - I*ky],[-I*kx + I*ky, 0]]),
                Matrix([[kz,  0],[ 0, kz]]),
                Matrix([[ 0, kz],[kz,  0]])
            ]
        )
    ),
])
def test_symmetric_hamiltonian(symmetry_operations, expr_basis, repr_basis, result):
    assert kp.symmetric_hamiltonian(
        *symmetry_operations, 
        expr_basis=expr_basis, 
        repr_basis=repr_basis
    ) == result

@pytest.mark.parametrize('symmetry_operations,expr_basis,repr_basis', [
    (
        [kp.SymmetryOperation(
            kmatrix=[[0, 1, 0], [1, 0, 0], [0, 0, 1]],
            repr=kp.Representation(
                matrix=[[0, 1], [1, 0]],
                complex_conjugate=False
            )
        )],
        kp.monomial_basis(0),
        'auto'
    ),
    (
        [kp.SymmetryOperation(
            kmatrix=[[0, 1, 0], [1, 0, 0], [0, 0, 1]],
            repr=kp.Representation(
                matrix=[[0, 1], [1, 0]],
                complex_conjugate=False
            )
        )],
        kp.monomial_basis(1),
        'auto'
    ),
    (
        [kp.SymmetryOperation(
            kmatrix=[[0, 1, 0], [1, 0, 0], [0, 0, 1]],
            repr=kp.Representation(
                matrix=[[0, 1], [1, 0]],
                complex_conjugate=False
            )
        )],
        kp.monomial_basis(*range(3)),
        'auto'
    ),
])
def test_symmetric_hamiltonian_consistency(symmetry_operations, expr_basis, repr_basis, compare_equal):
    compare_equal(sp.lambdify(
        (), 
        kp.symmetric_hamiltonian(
            *symmetry_operations, 
            expr_basis=expr_basis, 
            repr_basis=repr_basis
        )[0]
    )())
