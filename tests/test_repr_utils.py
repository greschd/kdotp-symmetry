#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import sympy as sp

from kdotp_symmetry._repr_utils import frobenius_product, hermitian_basis, hermitian_to_vector

sigma0 = sp.Matrix([[1, 0], [0, 1]])
sigmax = sp.Matrix([[0, 1], [1, 0]])
sigmay = sp.Matrix([[0, -sp.I], [sp.I, 0]])
sigmaz = sp.Matrix([[1, 0], [0, -1]])

sigma_vec = [sigma0, sigmax, sigmay, sigmaz]


@pytest.mark.parametrize('A', sigma_vec)
@pytest.mark.parametrize('B', sigma_vec)
def test_frobenius_product(A, B):
    if A == B:
        assert frobenius_product(A, B) == 2
    else:
        assert frobenius_product(A, B) == 0


@pytest.mark.parametrize(
    'dim,result', [(0, []), (1, [sp.Matrix([[1]])]), (
        2, [
            sp.Matrix([[1, 0], [0, 0]]),
            sp.Matrix([[0, 0], [0, 1]]),
            sp.Matrix([[0, 1], [1, 0]]),
            sp.Matrix([[0, -sp.I], [sp.I, 0]])
        ]
    ), (
        3, [
            sp.Matrix([[1, 0, 0], [0, 0, 0], [0, 0, 0]]),
            sp.Matrix([[0, 0, 0], [0, 1, 0], [0, 0, 0]]),
            sp.Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 1]]),
            sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
            sp.Matrix([[0, -sp.I, 0], [sp.I, 0, 0], [0, 0, 0]]),
            sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
            sp.Matrix([[0, 0, -sp.I], [0, 0, 0], [sp.I, 0, 0]]),
            sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
            sp.Matrix([[0, 0, 0], [0, 0, -sp.I], [0, sp.I, 0]]),
        ]
    )]
)
def test_hermitian_basis(dim, result):
    for b, r in zip(hermitian_basis(dim), result):
        print(b)
        print(r)
        assert b == r
    assert hermitian_basis(dim) == result


@pytest.mark.parametrize(
    'mat,vec,basis', [(
        sp.Matrix([[0, 1 + sp.I], [1 - sp.I, 0]]), (0, 0, 1, -1),
        hermitian_basis(2)
    ), (
        sp.Matrix([[2, sp.sqrt(2) + sp.I], [sp.sqrt(2) - sp.I, -3]]),
        (2, -3, sp.sqrt(2), -1), hermitian_basis(2)
    )]
)
def test_hermitian_to_vector(mat, vec, basis):
    assert hermitian_to_vector(mat, basis) == vec
