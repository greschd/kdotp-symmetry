# © 2017-2018, ETH Zurich, Institut für Theoretische Physik
# Author:  Dominik Gresch <greschd@gmx.ch>
"""
Defines the functions to calculate the basis of an intersection of vector spaces.
"""

import enum
from collections import namedtuple

import sympy as sp
import numpy as np
import networkx as nx

ZassenhausResult = namedtuple('ZassenhausResult', ['sum', 'intersection'])


def zassenhaus(basis_a, basis_b):
    r"""
    Given two bases ``basis_a`` and ``basis_b`` of vector spaces :math:`U, W \subseteq V`, computes bases of :math:`U + W` and :math:`U \cap W` using the Zassenhaus algorithm.
    """
    # handle the case where one of the bases is empty
    if len(basis_a) == 0 or len(basis_b) == 0:  # pylint: disable=no-else-return
        mat, pivot = sp.Matrix(basis_a).rref()
        plus_basis = mat[:len(pivot), :].tolist()
        return ZassenhausResult(sum=plus_basis, intersection=[])
    else:
        # Set up the Zassenhaus matrix from the given bases.
        A = sp.Matrix(basis_a)
        B = sp.Matrix(basis_b)
        dim = A.shape[1]  # pylint: disable=unsubscriptable-object
        if B.shape[1] != dim:  # pylint: disable=unsubscriptable-object
            raise ValueError('Inconsistent dimensions of the two bases given.')
        zassenhaus_mat = A.row_join(A).col_join(B.row_join(sp.zeros(*B.shape)))  # pylint: disable=not-an-iterable
        mat, pivot = zassenhaus_mat.rref()

        # idx is the row index of the first row belonging to the intersection basis
        idx = np.searchsorted(pivot, dim)
        plus_basis = mat[:idx, :dim].tolist()
        # The length of the pivot table is used to get rid of all-zero rows
        int_basis = mat[idx:len(pivot), dim:].tolist()
        return ZassenhausResult(sum=plus_basis, intersection=int_basis)


def intersection_basis(*bases):
    r"""
    Given ``bases`` of different subspaces :math:`U_i \subseteq V`, returns a basis of the intersection :math:`\bigcap_i U_i`. The basis vectors must all have the same length.
    """
    # sorting is not strictly needed, but usually more efficient
    bases_sorted = sorted(bases, key=len)
    current_basis = bases_sorted.pop(0)
    for basis in bases_sorted:
        if len(current_basis) == 0:
            break
        current_basis = zassenhaus(current_basis, basis).intersection
    return current_basis


@enum.unique
class _NodeType(enum.Enum):
    ROW = 'row'
    COLUMN = 'column'


def nullspace_blocked(matrix, **kwargs):
    """
    Calculate the nullspace of a given matrix. This is functionally equivalent to sympy's ``nullspace`` method, but it first subdivides the matrix into block-diagonal parts if possible.

    Keyword arguments are forwarded to the sympy nullspace method.
    """
    n_rows, n_cols = matrix.shape
    row_nodes = [(_NodeType.ROW, i) for i in range(n_rows)]
    column_nodes = [(_NodeType.COLUMN, i) for i in range(n_cols)]

    graph = nx.Graph()

    graph.add_nodes_from(row_nodes)
    graph.add_nodes_from(column_nodes)

    for row_idx in range(n_rows):
        for column_idx, val in enumerate(matrix.row(row_idx)):
            if val.is_nonzero:
                graph.add_edge((_NodeType.ROW, row_idx),
                               (_NodeType.COLUMN, column_idx))

    nullspace = []
    components = list(nx.connected_components(graph))
    for component in components:
        row_indices = sorted(
            idx for kind, idx in component if kind is _NodeType.ROW
        )
        column_indices = sorted(
            idx for kind, idx in component if kind is _NodeType.COLUMN
        )
        if len(column_indices) == 0:
            continue
        mat_part = matrix[row_indices, column_indices]
        # Get rid of fractions -- least common multiple of the denominators
        # This greatly improves the performance of sympy's nullspace -- for
        # whatever reason.
        mat_part *= sp.lcm([sp.fraction(val)[1] for val in mat_part])
        nullspace_part = np.array(mat_part.nullspace(**kwargs))
        if len(nullspace_part) == 0:
            continue
        nullspace_part_extended = np.zeros((len(nullspace_part), n_cols),
                                           dtype=object)
        nullspace_part_extended[:, column_indices] = nullspace_part
        nullspace.extend([sp.Matrix(vec) for vec in nullspace_part_extended])

    return nullspace
