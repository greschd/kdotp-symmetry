#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import sympy as sp

from kdotp_symmetry._expr_utils import expr_to_vector, monomial_basis, matrix_to_expr_operator

kx, ky, kz = sp.symbols('kx, ky, kz')


@pytest.mark.parametrize(
    'expr,vector,basis',
    [(1 + kx - ky + 2 * kz, (1, 1, -1, 2), [sp.Integer(1), kx, ky, kz]), (
        kx * ky + kx * ky * kz, (0, 0, 0, 0, 1, 0, 0, 1),
        [sp.Integer(1), kx, ky, kz, kx * ky, kx * kz, ky * kz, kx * ky * kz]
    ), (
        1 + sp.Rational(1, 4) * ky, (1, 0, 0.25, 0),
        [sp.Integer(1), kx, ky, kz]
    )]
)
def test_expr_to_vector(expr, vector, basis):
    assert expr_to_vector(expr, basis=basis) == vector


@pytest.mark.parametrize(
    'expr,basis', [
        (1 + kx, [sp.Integer(1), kx, kx, kz]),
    ]
)
def test_basis_not_independent(expr, basis):
    with pytest.raises(ValueError):
        expr_to_vector(expr, basis=basis)


@pytest.mark.parametrize(
    'dim,basis', [(0, [sp.Integer(1)]), (1, [sp.Integer(1), kx, ky, kz]), (
        2, [
            sp.Integer(1), kx, ky, kz, kx**2, kx * ky, kx * kz, ky**2, ky * kz,
            kz**2
        ]
    )]
)
def test_monomial_basis(dim, basis):
    assert monomial_basis(*range(dim + 1)) == basis


def test_monomial_basis_negative_degree():
    with pytest.raises(ValueError):
        monomial_basis(1, 2, -3)


@pytest.mark.parametrize(
    'matrix_form,expr1,expr2',
    [([[0, 1, 0], [1, 0, 0], [0, 0, -1]], 2 + kx**2 + kx * ky + kx * kz,
      2 + ky**2 + kx * ky - ky * kz),
     ([[0, 1, 0], [1, 0, 0], [0, 0, -1]], -1 + kx**2 * ky + kx - kx * kz**2,
      -1 + ky**2 * kx + ky - ky * kz**2)]
)
def test_matrix_to_expr_operator(matrix_form, expr1, expr2):
    assert sp.simplify(
        sp.Eq(matrix_to_expr_operator(matrix_form)(expr1), expr2)
    )


@pytest.mark.parametrize(
    'matrix_form,expr1,expr2',
    [([[0, 1, 0], [1, 0, 0], [0, 0, -1]], 2 + kx**2 + kx * ky + kx * kz,
      2 + ky**2 + kx * ky - ky * kz),
     ([[0, 1, 0], [1, 0, 0], [0, 0, -1]], -1 + kx**2 * ky + kx - kx * kz**2,
      -1 + ky**2 * kx + ky - ky * kz**2), ([[0, 1, 0], [0, 0, 1], [1, 0, 0]],
                                           2 + kx**2 + kx * ky + kx * kz,
                                           2 + kz**2 + kz * kx + kz * ky)]
)
def test_matrix_to_expr_operator_double_eval(matrix_form, expr1, expr2):
    op = matrix_to_expr_operator(matrix_form)
    assert sp.simplify(sp.Eq(op(expr1), expr2))
    assert sp.simplify(sp.Eq(op(expr1), expr2))
