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
from fsc.export import export

K_VEC = sp.symbols('kx, ky, kz')

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
    vec = tuple(v.nsimplify() for v in vec)
    # check consistency
    assert expr.equals(sum(v * b for v, b in zip(vec, basis)))
    return vec

@export
def monomial_basis(*degrees):
    """
    Returns the product basis of (kx, ky, kz), with monomials of the given degrees.

    :param degrees: Degree of the monomials. Multiple degrees can be given, in which case the basis consists of the monomials of all given degrees.
    :type degrees: int

    Example:

        >>> import kdotp_symmetry as kp
        >>> kp.monomial_basis(*range(3))
        [1, kx, ky, kz, kx**2, kx*ky, kx*kz, ky**2, ky*kz, kz**2]
    """
    if any(p < 0 for p in degrees):
        raise ValueError('Degrees must be non-negative integers')
    basis = []
    for d in sorted(degrees):
        monomial_tuples = combinations_with_replacement(K_VEC, d)
        basis.extend(
            reduce(operator.mul, m, sp.Integer(1))
            for m in monomial_tuples
        )
    return basis

def matrix_to_expr_operator(matrix_form, repr_has_cc=False):
    """Returns a function that operates on expression, corresponding to the given ``matrix_form`` which operates on a vector in real space. ``repr_has_cc`` determines whether the symmetry contains time reversal."""
    # k-form and r-form of the matrix are related by A -> A^-1^T
    # => matrix^T gives g^-1 in k-space coordinates
    # Change sign if the representation has complex conjugation
    k_matrix_form = sp.Matrix(matrix_form).T
    if repr_has_cc:
        k_matrix_form *= -1
    substitution = list(zip(K_VEC, k_matrix_form @ sp.Matrix(K_VEC)))
    def operator(expr):
        return expr.subs(substitution, simultaneous=True)
    return operator
