#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    21.01.2017 09:41:18 CET
# File:    expr_to_vector.py

import random
import operator
from functools import reduce
from itertools import combinations_with_replacement

import sympy as sp

from .constants import K_VEC

def expr_to_vector(
        expr,
        basis,
        *,
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
        vals = [(k, random_fct()) for k in K_VEC]
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
    assert sp.simplify(sp.Eq(sum(v * b for v, b in zip(vec, basis)), expr))
    return vec

def monomial_basis(power):
    """Returns the monomial basis of (kx, ky, kz), up to the chosen ``power``."""
    if power < 0:
        raise ValueError('The power must be a non-negative integer.')
    basis = []
    for p in range(power + 1):
        monomial_tuples = combinations_with_replacement(K_VEC, p)
        basis.extend(
            reduce(operator.mul, m, sp.Integer(1))
            for m in monomial_tuples
        )
    return basis
    
def matrix_to_expr_operator(k_matrix_form):
    """Returns a function that operates on expression, corresponding to the given ``k_matrix_form`` which operates on a vector in k-space."""
    substitution = list(zip(
        K_VEC, 
        list(
            sp.Matrix(k_matrix_form) @ sp.Matrix(K_VEC)
        )
    ))
    def operator(expr):
        return expr.subs(substitution, simultaneous=True)
    return operator
