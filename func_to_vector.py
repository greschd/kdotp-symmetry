#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>

import random
import sympy as sp

def func_to_vector(
        expr,
        basis,
        *,
        coordinates=sp.symbols('kx, ky, kz'),
        random_fct=lambda: random.randint(-100, 100)
    ):
    dim = len(basis)
    A = []
    b = []
    for _ in range(2 * dim):
        if sp.Matrix(A).rank() >= len(basis):
            break
        vals = [(k, random_fct()) for k in coordinates]
        A.append([b.subs(vals) for b in basis])
        b.append(expr.subs(vals))
    else:
        raise ValueError('Could not find a sufficient number of linearly independent vectors')

    res = sp.linsolve((sp.Matrix(A), sp.Matrix(b)), sp.symbols('a b c'))
    if len(res) != 1:
        raise ValueError('No or multiple results found: {}'.format(res))
    return next(iter(res))

if __name__ == '__main__':
    kx, ky, kz = sp.symbols('kx, ky, kz')

    expr = kx * (4 + ky) + 2 * ky + 3 # [3, 4, 2, 1]
    basis = [sp.Integer(1), kx, ky, kx * ky]
    print(to_vector(expr, basis))
