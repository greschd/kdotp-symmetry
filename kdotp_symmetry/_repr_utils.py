# © 2017-2018, ETH Zurich, Institut für Theoretische Physik
# Author:  Dominik Gresch <greschd@gmx.ch>
"""
Utilities for handling symmetry representations, such as converting them to matrix form.
"""

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
        basis.append(sp.SparseMatrix(dim, dim, {(i, i): 1}))

    # off-diagonal entries
    for i in range(dim):
        for j in range(i + 1, dim):
            # real
            basis.append(sp.SparseMatrix(dim, dim, {(i, j): 1, (j, i): 1}))

            # imag
            basis.append(
                sp.SparseMatrix(dim, dim, {(i, j): -sp.I,
                                           (j, i): sp.I})
            )

    assert len(basis) == dim**2

    return basis


def check_orthogonal(basis):
    """Check orthogonality for a given ``basis``."""
    for i, b_i in enumerate(basis):
        for offset, b_j in enumerate(basis[i:]):
            j = i + offset
            frob_product = frobenius_product(b_i, b_j)
            if i == j:
                if frob_product != 1:
                    raise ValueError(
                        'Basis element {} has norm {}, not one.'.format(
                            i, frob_product
                        )
                    )
            else:
                if frob_product != 0:
                    raise ValueError(
                        'Basis elements {}, {} are not orthogonal.'.format(
                            i, j
                        )
                    )


def hermitian_to_vector(matrix, basis, basis_norm_squares=None):
    """
    Returns a the vector representing the ``matrix`` w.r.t. the given *orthogonal* ``basis``.
    """
    vec = tuple(
        frobenius_product(matrix, b) / norm_sq for b, norm_sq in zip(
            basis, basis_norm_squares
            or [frobenius_product(b, b) for b in basis]
        )
    )
    vec = tuple(v.nsimplify() for v in vec)
    # check consistency
    if not matrix.equals(
        sum((v * b for v, b in zip(vec, basis)), sp.zeros(*matrix.shape))
    ):
        raise ValueError(
            'Vector {vec} in basis {basis} does not match matrix {matrix}'.
            format(vec=vec, basis=basis, matrix=matrix)
        )
    return vec


def repr_to_matrix_operator(matrix_representation, complex_conjugate=False):
    """
    Converts a symmetry representation into the corresponding matrix operator.

    :param matrix_representation: Real-space matrix form of the symmetry representation.
    :type matrix_representation: sympy.Matrix

    :param complex_conjugate: Specifies whether the representation contains complex conjugation.
    :type complex_conjugate: bool
    """
    matrix_representation = sp.Matrix(matrix_representation)

    def operator(matrix):
        if complex_conjugate:
            return matrix_representation @ matrix.conjugate(
            ) @ matrix_representation.H
        return matrix_representation @ matrix @ matrix_representation.H

    return operator
