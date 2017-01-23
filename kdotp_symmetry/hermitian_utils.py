# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    21.01.2017 09:41:11 CET
# File:    matrix_to_vector.py

import sympy as sp

def frobenius_product(A, B):
    r"""
    Returns the Frobenius scalar product <A, B> = Tr(A^\dagger B) for two matrices.
    """
    return (A.conjugate().transpose() @ B).trace().simplify()

def create_onb_hermitian(dim):
    """
    Returns an orthonormal basis (w.r.t. the Frobenius scalar product) of the hermitian matrices of size 'dim'.
    """
    basis = []
    # diagonal entries
    for i in range(dim):
        mat = sp.zeros(dim)
        mat[i, i] = 1
        basis.append(mat)

    # off-diagonal entries
    x = 1 / sp.sqrt(2)
    for i in range(dim):
        for j in range(i + 1, dim):
            # real
            mat = sp.zeros(dim)
            mat[i, j] = x
            mat[j, i] = x
            basis.append(mat)

            # imag
            mat = sp.zeros(dim)
            mat[i, j] = -1j * x
            mat[j, i] = 1j * x
            basis.append(mat)

    # check ONB property
    assert len(basis) == dim**2
    _assert_onb(basis)

    return basis

def _assert_onb(basis):
    """Check ONB properties for a given ``basis``."""
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            if i == j:
                assert frobenius_product(bi, bj) == 1
            else:
                assert frobenius_product(bi, bj) == 0

def hermitian_to_vector(matrix, basis):
    """
    Returns a the vector representing the ``matrix`` w.r.t. the given *orthonormal* ``basis``.
    """
    _assert_onb(basis)
    vec = tuple(frobenius_product(matrix, b) for b in basis)
    # check consistency
    assert sp.Eq(sum((v * b for v, b in zip(vec, basis)), sp.zeros(*matrix.shape)), matrix)
    return vec
