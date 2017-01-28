#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    25.01.2017 15:14:49 CET
# File:    analyse.py

from collections import namedtuple

import sympy as sp
from sympy.physics.quantum import TensorProduct
import numpy as np
from fsc.export import export

from ._expr_utils import expr_to_vector, monomial_basis, matrix_to_expr_operator
from ._repr_utils import hermitian_to_vector, hermitian_basis, repr_to_matrix_operator
from ._linalg import intersection_basis
from ._to_matrix import to_matrix

SymmetryOperation = namedtuple('SymmetryOperation', ['kmatrix', 'repr'])
SymmetryOperation.__doc__ = 'TODO'
SymmetryOperation.kmatrix.__doc__ = 'TODO'
SymmetryOperation.repr.__doc__ = 'TODO'
Representation = namedtuple('Representation', ['matrix', 'complex_conjugate'])
Representation.__doc__ = 'TODO'
Representation.matrix.__doc__ = 'TODO'
Representation.complex_conjugate.__doc__ = 'TODO'

__all__ = ['Representation', 'SymmetryOperation']

@export
def symmetric_hamiltonian(*symmetry_operations, expr_basis, repr_basis='auto'):
    """..."""
    expr_dim = len(expr_basis)
    # for sympy or numpy matrices
    try:
        repr_matrix_size = symmetry_operations[0].repr.matrix.shape[0]
    # for plain lists -- this doesn't work for sympy matrices because
    # their 'len' is the total number of elements
    except AttributeError:
        repr_matrix_size = len(symmetry_operations[0].repr.matrix)

    if repr_basis == 'auto':
        repr_basis = hermitian_basis(repr_matrix_size)
    repr_dim = len(repr_basis)
    full_dim = expr_dim * repr_dim
    full_basis = [
        sp.Matrix(x) for x in 
        np.outer(expr_basis, repr_basis).reshape(full_dim, repr_matrix_size, repr_matrix_size).tolist()
    ]
    
    invariant_bases = []
    for op in symmetry_operations:
        # create the matrix form of the two operators
        expr_mat = to_matrix(
            operator=matrix_to_expr_operator(op.kmatrix),
            basis=expr_basis,
            to_vector_fct=expr_to_vector
        )
        repr_mat = to_matrix(
            operator=repr_to_matrix_operator(*op.repr),
            basis=repr_basis,
            to_vector_fct=hermitian_to_vector
        )
        # outer product
        full_mat = TensorProduct(expr_mat, repr_mat)

        # get Eig(F \ocross G, 1) basis
        invariant_bases.append(
            np.array(
                (full_mat - sp.eye(full_dim)).nullspace()
            ).tolist()
        )
    
    basis_vectors = intersection_basis(*invariant_bases)
    basis_vectors_expanded = []
    for vec in basis_vectors:
        basis_vectors_expanded.append(
            sum((v * b for v, b in zip(vec, full_basis)), sp.zeros(repr_matrix_size))
        )
    return basis_vectors_expanded

