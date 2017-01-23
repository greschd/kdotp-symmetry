#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    23.01.2017 10:55:46 CET
# File:    test_func_utils.py

from collections import namedtuple

import pytest
import sympy as sp

from kdotp_symmetry.func_utils import func_to_vector, create_monomial_basis, operator_form

kx, ky, kz = sp.symbols('kx, ky, kz')
Expression = namedtuple('Expression', ('expr', 'vector', 'basis'))

expressions = [
    Expression(
        1 + kx - ky + 2 * kz,
        (1, 1, -1, 2),
        [sp.Integer(1), kx, ky, kz]
    ),
    Expression(
        kx * ky + kx * ky * kz,
        (0, 0, 0, 0, 1, 0, 0, 1),
        [sp.Integer(1), kx, ky, kz, kx * ky, kx * kz, ky * kz, kx * ky * kz]
    ),
    Expression(
        1 + 0.25 * ky,
        (1, 0, 0.25, 0),
        [sp.Integer(1), kx, ky, kz]
    )
]

@pytest.mark.parametrize('expr', expressions)
def test_func_to_vector(expr):
    assert func_to_vector(expr.expr, basis=expr.basis) == expr.vector
