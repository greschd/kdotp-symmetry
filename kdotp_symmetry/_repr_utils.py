# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    21.01.2017 09:41:11 CET
# File:    matrix_to_vector.py

import sympy as sp
from fsc.export import export

def frobenius_product(A, B):
    r"""
    Returns the Frobenius scalar product <A, B> = Tr(A^\dagger B) for two matrices.
    """
    return (A.H @ B).trace().simplify()

@export
def hermitian_basis(dim):
    """
    Returns a basis of the hermitian matrices of size ``dim`` that is orthogonal w.r.t. the Frobenius scalar product.

    :param dim: size of the matrices
    :type dim:  int

    Example:

        >>> import kdotp_symmetry as kp
        >>> kp.hermitian_basis(2)
        [Matrix([
        [1, 0],
        [0, 0]]), Matrix([
        [0, 0],
        [0, 1]]), Matrix([
        [0, 1],
        [1, 0]]), Matrix([
        [0, -I],
        [I,  0]])]
    """
    basis = []
    # diagonal entries
    for i in range(dim):
        mat = sp.zeros(dim)
        mat[i, i] = 1
        basis.append(mat)

    # off-diagonal entries
    for i in range(dim):
        for j in range(i + 1, dim):
            # real
            mat = sp.zeros(dim)
            mat[i, j] = 1
            mat[j, i] = 1
            basis.append(mat)

            # imag
            mat = sp.zeros(dim)
            mat[i, j] = -sp.numbers.I
            mat[j, i] = sp.numbers.I
            basis.append(mat)

    # check ONB property
    assert len(basis) == dim**2
    _assert_orthogonal(basis)

    return basis

def _assert_orthogonal(basis):
    """Check orthogonality for a given ``basis``."""
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            if i == j:
                assert frobenius_product(bi, bj) != 0
            else:
                assert frobenius_product(bi, bj) == 0

def hermitian_to_vector(matrix, basis):
    """
    Returns a the vector representing the ``matrix`` w.r.t. the given *orthonormal* ``basis``.
    """
    _assert_orthogonal(basis)
    vec = tuple(frobenius_product(matrix, b) / frobenius_product(b, b) for b in basis)
    vec = tuple(v.nsimplify() for v in vec)
    # check consistency
    assert matrix.equals(sum((v * b for v, b in zip(vec, basis)), sp.zeros(*matrix.shape)))
    return vec

def repr_to_matrix_operator(matrix_representation, complex_conjugate=False):
    matrix_representation = sp.Matrix(matrix_representation)
    def operator(matrix):
        if complex_conjugate:
            return matrix_representation @  matrix.conjugate() @ matrix_representation.H
        else:
            return matrix_representation @  matrix @ matrix_representation.H
    return operator
