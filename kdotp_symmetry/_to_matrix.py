# © 2017-2018, ETH Zurich, Institut für Theoretische Physik
# Author:  Dominik Gresch <greschd@gmx.ch>
"""
Defines a function to convert an operator into matrix form.
"""

import types

import sympy as sp


def to_matrix(
    *,
    operator,
    basis,
    to_vector_fct,
    to_vector_kwargs=types.MappingProxyType({})
):
    """
    Convert an operator into matrix form, w.r.t. a given basis.

    :param operator: Operator which is to be expressed in matrix form.

    :param basis: Basis with respect to which the matrix should be written.
    :type basis: list

    :param to_vector_fct: Function which turns elements in the vector space into vector form.

    :param to_vector_kwargs: Additional keyword arguments passed to ``to_vector_fct``.
    :type to_vector_kwargs: Mapping
    """
    return sp.Matrix([
        to_vector_fct(operator(b), basis=basis, **to_vector_kwargs)
        for b in basis
    ]).transpose()
