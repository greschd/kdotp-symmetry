#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    23.01.2017 14:41:35 CET
# File:    test_matrix_representation.py

import pytest
import sympy as sp

from kdotp_symmetry.matrix_representation import matrix_representation

from kdotp_symmetry.constants import K_VEC
from kdotp_symmetry.func_utils import func_to_vector, create_monomial_basis, operator_form

@pytest.mark.parametrize('operator,basis,to_vector_fct,result', [
    (
        operator_form([[0, 1, 0], [1, 0, 0], [0, 0, -1]]),
        create_monomial_basis(2),
        func_to_vector,
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
    )
])
def test_matrix_representation(operator, basis, to_vector_fct, result):
    kx, ky, kz = K_VEC
    #~ print([operator(b) for b in basis])
    print(operator(kx * kz))
    print(matrix_representation(operator, basis, to_vector_fct))
    assert matrix_representation(operator, basis, to_vector_fct) == result
