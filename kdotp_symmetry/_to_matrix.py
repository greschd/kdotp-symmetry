"""
Defines a function to convert an operator into matrix form.
"""

import sympy as sp


def to_matrix(operator, basis, to_vector_fct):
    """
    Convert an operator into matrix form, w.r.t. a given basis.

    :param operator: Operator which is to be expressed in matrix form.

    :param basis: Basis with respect to which the matrix should be written.
    :type basis: list

    :param to_vector_fct: Function which turns elements in the vector space into vector form.
    """
    return sp.Matrix([to_vector_fct(operator(b), basis=basis)
                      for b in basis]).transpose()
