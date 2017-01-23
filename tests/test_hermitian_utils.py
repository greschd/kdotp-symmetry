#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    23.01.2017 13:56:42 CET
# File:    test_hermitian_utils.py

import pytest
import sympy as sp

from kdotp_symmetry.hermitian_utils import frobenius_product, create_onb_hermitian, hermitian_to_vector

sigma0 = sp.Matrix([[1, 0], [0, 1]])
sigmax = sp.Matrix([[0, 1], [1, 0]])
sigmay = sp.Matrix([[0, -1j], [1j, 0]])
sigmaz = sp.Matrix([[1, 0], [0, -1]])

sigma_vec = [sigma0, sigmax, sigmay, sigmaz]

@pytest.mark.parametrize('A', sigma_vec)
@pytest.mark.parametrize('B', sigma_vec)
def test_frobenius_product(A, B):
    if A == B:
        assert frobenius_product(A, B) == 2
    else:
        assert frobenius_product(A, B) == 0

@pytest.mark.parametrize('dim,result', [
    (0, []),
    (1, [sp.Matrix([[1]])]),
    (2, [
        sp.Matrix([[1, 0], [0, 0]]), 
        sp.Matrix([[0, 0], [0, 1]]),
        sp.Integer(1) / sp.sqrt(2) * sp.Matrix([[0, 1], [1, 0]]),
        sp.Integer(1) / sp.sqrt(2) * sp.Matrix([[0, -1j], [1j, 0]])
    ]),
    (3, [
        sp.Matrix([[1, 0, 0], [0, 0, 0], [0, 0, 0]]),
        sp.Matrix([[0, 0, 0], [0, 1, 0], [0, 0, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 1]]),
        1 / sp.sqrt(2) * sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
        1 / sp.sqrt(2) * sp.Matrix([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]]),
        1 / sp.sqrt(2) * sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
        1 / sp.sqrt(2) * sp.Matrix([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]]),
        1 / sp.sqrt(2) * sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
        1 / sp.sqrt(2) * sp.Matrix([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]]),
    ])
])
def test_create_onb_hermitian(dim, result):
    for b, r in zip(create_onb_hermitian(dim), result):
        print(b)
        print(r)
        assert b == r
    assert create_onb_hermitian(dim) == result

@pytest.mark.parametrize('mat,vec,basis', [
    (
        sp.Matrix([[0, 1 + 1j], [1 - 1j, 0]]), 
        (0, 0, sp.sqrt(2), -sp.sqrt(2)), 
        create_onb_hermitian(2)
    ),
    (
        sp.Matrix([[2, sp.sqrt(2) + 1j], [sp.sqrt(2) - 1j, -3]]), 
        (2, -3, 2, -sp.sqrt(2)), 
        create_onb_hermitian(2)
    )
])
def test_hermitian_to_vector(mat, vec, basis):
    assert hermitian_to_vector(mat, basis) == vec
