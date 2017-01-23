#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    23.01.2017 10:55:46 CET
# File:    test_func_utils.py

import pytest
import sympy as sp

from kdotp_symmetry.func_utils import func_to_vector, create_monomial_basis, operator_form

kx, ky, kz = sp.symbols('kx, ky, kz')

@pytest.mark.parametrize('expr,vector,basis', [
    (
        1 + kx - ky + 2 * kz,
        (1, 1, -1, 2),
        [sp.Integer(1), kx, ky, kz]
    ),
    (
        kx * ky + kx * ky * kz,
        (0, 0, 0, 0, 1, 0, 0, 1),
        [sp.Integer(1), kx, ky, kz, kx * ky, kx * kz, ky * kz, kx * ky * kz]
    ),
    (
        1 + sp.Rational(1, 4) * ky,
        (1, 0, 0.25, 0),
        [sp.Integer(1), kx, ky, kz]
    )
])
def test_func_to_vector(expr, vector, basis):
    assert func_to_vector(expr, basis=basis) == vector

@pytest.mark.parametrize('dim,basis', [
    (0, [sp.Integer(1)]),
    (1, [sp.Integer(1), kx, ky, kz]),
    (2, [sp.Integer(1), kx, ky, kz, kx**2, kx * ky, kx * kz, ky**2, ky * kz, kz**2])
])
def test_monomial_basis(dim, basis):
    assert create_monomial_basis(dim) == basis
