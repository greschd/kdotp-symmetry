# © 2017-2018, ETH Zurich, Institut für Theoretische Physik
# Author:  Dominik Gresch <greschd@gmx.ch>
"""
Defines functions to construct the basis of the symmetry-constrained Hamiltonian.
"""

import sympy as sp
from sympy.physics.quantum import TensorProduct
import numpy as np
import scipy.linalg as la
from fsc.export import export

from ._expr_utils import expr_to_vector, matrix_to_expr_operator
from ._repr_utils import hermitian_to_vector, hermitian_basis, repr_to_matrix_operator, check_orthogonal, frobenius_product
from ._linalg import intersection_basis, nullspace_blocked
from ._to_matrix import to_matrix
from ._logging_setup import LOGGER


@export
def symmetric_hamiltonian(
    *symmetry_operations,
    expr_basis,
    repr_basis='auto',
    check_repr_basis=False
):
    r"""
    Calculates the basis of the symmetric Hamiltonian for a given set of symmetry operations.

    :param symmetry_operations: The symmetry operations that the Hamiltonian should respect.
    :type symmetry_operations: :py:class:`symmetry_representation.SymmetryOperation`

    :param expr_basis: The basis for the :math:`\mathbf{k}`-functions that are considered.
    :type expr_basis: :py:class:`list` of :py:mod:`sympy` expressions

    :param repr_basis: The basis for the hermitian matrices, with the same size as the representations. By default, the :py:func:`.hermitian_basis` of the appropriate size is used.
    :type repr_basis: :py:class:`list` of :py:mod:`sympy` matrices

    :param check_repr_basis: Flag to enable explicitly checking the orthogonality of ``repr_basis``.
    :type check_repr_basis: bool

    :returns: Basis for the symmetric Hamiltonian, as a :py:class:`list` of :py:mod:`sympy` matrix expressions.
    """
    if any(sym_op.numeric for sym_op in symmetry_operations):
        raise ValueError(
            'Symmetry operations used in kdotp-symmetry can not be numeric.'
        )
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

    if check_repr_basis:
        check_orthogonal(repr_basis)

    repr_dim = len(repr_basis)
    repr_basis_norm_squares = [frobenius_product(b, b) for b in repr_basis]
    full_dim = expr_dim * repr_dim
    full_basis = [
        sp.Matrix(x) for x in np.outer(expr_basis, repr_basis).
        reshape(full_dim, repr_matrix_size, repr_matrix_size).tolist()
    ]

    invariant_bases = []
    for sym_op in symmetry_operations:
        LOGGER.info('Calculating matrix form of expression.')
        expr_mat = to_matrix(
            operator=matrix_to_expr_operator(
                sym_op.rotation_matrix, repr_has_cc=sym_op.repr.has_cc
            ),
            basis=expr_basis,
            to_vector_fct=expr_to_vector
        )
        LOGGER.info('Calculating matrix form of representation.')
        repr_mat = to_matrix(
            operator=repr_to_matrix_operator(
                sym_op.repr.matrix, complex_conjugate=sym_op.repr.has_cc
            ),
            basis=repr_basis,
            to_vector_fct=hermitian_to_vector,
            to_vector_kwargs=dict(basis_norm_squares=repr_basis_norm_squares)
        )
        # outer product
        LOGGER.info('Calculating outer product.')
        full_mat = TensorProduct(expr_mat, repr_mat)

        # get Eig(F \ocross G, 1) basis
        mat = full_mat - sp.eye(full_dim)
        LOGGER.info('Calculating nullspace.')
        curr_basis = np.array(nullspace_blocked(mat, simplify=sp.nsimplify)
                              ).tolist()
        if len(curr_basis) != _numeric_nullspace_dim(mat):
            raise ValueError(
                'Analytic and numeric dimensions of the nullspace of the matrix {mat} do not match'
                .format(mat=mat)
            )
        invariant_bases.append(curr_basis)

    LOGGER.info('Calculating basis intersection.')
    basis_vectors = intersection_basis(*invariant_bases)
    LOGGER.info('Expanding basis vectors.')
    basis_vectors_expanded = []
    for vec in basis_vectors:
        basis_vectors_expanded.append(
            sum((v * b for v, b in zip(vec, full_basis)),
                sp.zeros(repr_matrix_size))
        )
    return basis_vectors_expanded


def _numeric_nullspace_dim(mat):
    """Numerically computes the nullspace dimension of a matrix."""
    mat_numeric = np.array(mat.evalf().tolist(), dtype=complex)
    eigenvals = la.eigvals(mat_numeric)
    return np.sum(np.isclose(eigenvals, np.zeros_like(eigenvals)))
