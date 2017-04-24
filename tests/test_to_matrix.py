#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    23.01.2017 14:41:35 CET
# File:    test_to_matrix.py

import pytest
import sympy as sp

from kdotp_symmetry._to_matrix import to_matrix

from kdotp_symmetry._expr_utils import K_VEC, expr_to_vector, monomial_basis, matrix_to_expr_operator

from kdotp_symmetry._repr_utils import hermitian_to_vector, hermitian_basis, repr_to_matrix_operator

@pytest.mark.parametrize('operator,basis,to_vector_fct,result', [
    (
        matrix_to_expr_operator([[0, 1, 0], [1, 0, 0], [0, 0, -1]]),
        monomial_basis(*range(3)),
        expr_to_vector,
        sp.Matrix([
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, -1, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        ])
    ),
    (
        matrix_to_expr_operator([[0, 1, 0], [1, 0, 0], [0, 0, -1]]),
        monomial_basis(*range(2)),
        expr_to_vector,
        sp.Matrix([
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 1, 0, 0],
            [0, 0, 0, -1]
        ])
    ),
    (
        repr_to_matrix_operator(sp.Matrix([[0, 1], [1, 0]])),
        hermitian_basis(2),
        hermitian_to_vector,
        sp.Matrix([
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, -1],
        ])
    ),
    (
        repr_to_matrix_operator(sp.Matrix([[0, 1], [1, 0]]), complex_conjugate=True),
        hermitian_basis(2),
        hermitian_to_vector,
        sp.Matrix([
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ])
    ),
    (
        repr_to_matrix_operator(sp.eye(2), complex_conjugate=True),
        hermitian_basis(2),
        hermitian_to_vector,
        sp.Matrix([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, -1],
        ])
    ),
    (
        repr_to_matrix_operator(
            sp.diag(
                sp.exp(-sp.I * sp.pi * sp.Rational(2, 3)),
                sp.exp(sp.I * sp.pi * sp.Rational(2, 3))
            ),
            complex_conjugate=False
        ),
        hermitian_basis(2),
        hermitian_to_vector,
        sp.Matrix([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, -sp.Rational(1, 2), sp.sqrt(3) / 2],
            [0, 0, -sp.sqrt(3) / 2, -sp.Rational(1, 2)],
        ])
    )
])
def test_to_matrix(operator, basis, to_vector_fct, result):
    # kx, ky, kz = K_VEC
    print(to_matrix(operator, basis, to_vector_fct))
    assert to_matrix(operator, basis, to_vector_fct) == result
