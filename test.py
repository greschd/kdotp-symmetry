#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    23.01.2017 10:41:07 CET
# File:    test.py

from kdotp_symmetry.func_utils import *

if __name__ == '__main__':
    kx, ky, kz = sp.symbols('kx, ky, kz')
    expr = kx * (4 + ky) + 2 * ky + 3 # [3, 4, 2, 1]
    basis = [sp.Integer(1), kx, ky, kx * ky]
    print(func_to_vector(expr, basis))
    print(create_monomial_basis(3))
    
    op = operator_form([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
    print(op(kx * ky + ky + kz))
