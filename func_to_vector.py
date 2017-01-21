#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    21.01.2017 09:41:18 CET
# File:    func_to_vector.py

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
    # create random values for the coordinates and evaluate 
    # both the basis functions and the expression to generate
    # the linear equation to be solved
    A = []
    b = []
    for _ in range(2 * dim):
        if sp.Matrix(A).rank() >= len(basis):
            break
        vals = [(k, random_fct()) for k in coordinates]
        A.append([b.subs(vals) for b in basis])
        b.append(expr.subs(vals))
    else:
        # this could happen if the random_fct is bad, or the 'basis' is not
        # linearly independent
        raise ValueError('Could not find a sufficient number of linearly independent vectors')

    res = sp.linsolve((sp.Matrix(A), sp.Matrix(b)), sp.symbols('a b c'))
    if len(res) != 1:
        raise ValueError('No or multiple results found: {}'.format(res))
    vec = next(iter(res))
    # check consistency
    assert sp.Eq(sum(v * b for v, b in zip(vec, basis)), expr).simplify()
    return vec

if __name__ == '__main__':
    kx, ky, kz = sp.symbols('kx, ky, kz')

    expr = kx * (4 + ky) + 2 * ky + 3 # [3, 4, 2, 1]
    basis = [sp.Integer(1), kx, ky, kx * ky]
    print(func_to_vector(expr, basis))
