#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import sympy as sp
from sympy import Matrix
from sympy.core.numbers import I
import sympy.physics.matrices as sm
from sympy.physics.quantum import TensorProduct
import symmetry_representation as sr

import kdotp_symmetry as kp

kx, ky, kz = sp.symbols('kx, ky, kz')
PAULI_VEC = [sp.eye(2), *(sm.msigma(i) for i in range(1, 4))]

@pytest.mark.parametrize('symmetry_operations,expr_basis,repr_basis,result', [
    (
        [sr.SymmetryOperation(
            rotation_matrix=[[0, 1, 0], [1, 0, 0], [0, 0, 1]],
            repr_matrix=[[0, 1], [1, 0]],
            repr_has_cc=False
        )],
        kp.monomial_basis(0),
        'auto',
        [
            Matrix([[1, 0], [0, 1]]),
            Matrix([[0, 1], [1, 0]])
        ]
    ),
    (
        [sr.SymmetryOperation(
            rotation_matrix=[[0, 1, 0], [1, 0, 0], [0, 0, 1]],
            repr_matrix=[[0, 1], [1, 0]],
            repr_has_cc=False
        )],
        kp.monomial_basis(1),
        'auto',
        [
            Matrix([[ky, 0],[0, kx]]),
            Matrix([[kx,  0],[0, ky]]),
            Matrix([[0, kx + ky],[kx + ky, 0]]),
            Matrix([[0, I*kx - I*ky],[-I*kx + I*ky, 0]]),
            Matrix([[kz,  0],[ 0, kz]]),
            Matrix([[ 0, kz],[kz,  0]])
        ]
    ),
    (
        [
            sr.SymmetryOperation(
                rotation_matrix=[[0, 1, 0], [1, 0, 0], [0, 0, -1]],
                repr_matrix=sp.diag(I, -I, I, -I),
                repr_has_cc=False
            ),
            sr.SymmetryOperation(
                rotation_matrix=-sp.eye(3),
                repr_matrix=sp.diag(1, 1, -1, -1),
                repr_has_cc=False
            ),
            sr.SymmetryOperation(
                rotation_matrix=sp.eye(3),
                repr_matrix=TensorProduct(sp.eye(2), sp.Matrix([[0, -1], [1, 0]])),
                repr_has_cc=True
            )
        ],
        kp.monomial_basis(0),
        [TensorProduct(p1, p2) for p1 in PAULI_VEC for p2 in PAULI_VEC],
        [
            Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
            Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]),
        ]
    )
])
def test_symmetric_hamiltonian(symmetry_operations, expr_basis, repr_basis, result):
    assert kp.symmetric_hamiltonian(
        *symmetry_operations,
        expr_basis=expr_basis,
        repr_basis=repr_basis
    ) == result
