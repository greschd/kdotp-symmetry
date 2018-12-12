#!/usr/bin/env python
# -*- coding: utf-8 -*-

# © 2017-2018, ETH Zurich, Institut für Theoretische Physik
# Author:  Dominik Gresch <greschd@gmx.ch>

import sympy as sp
import symmetry_representation as sr

import kdotp_symmetry as kp

orbitals = [
    sr.Orbital(position=coord, function_string=fct, spin=spin)
    for spin in (sr.SPIN_UP, sr.SPIN_DOWN) for coord in
    (sp.Matrix([sp.Rational(1, 2)] * 3), sp.Matrix([sp.Rational(3, 4)] * 3))
    for fct in sr.WANNIER_ORBITALS['sp3']
]

symmetry_generators = [
    sr.SymmetryOperation.from_orbitals(
        orbitals=orbitals,
        real_space_operator=sr.RealSpaceOperator(
            rotation_matrix=sp.Matrix([[1, 1, 1], [0, -1, 0], [0, 0, -1]]),
            translation_vector=sp.Matrix([sp.Rational(1, 4)] * 3),
        ),
        rotation_matrix_cartesian=sp.Matrix([[-1, 0, 0], [0, 0, 1], [0, 1,
                                                                     0]]),
        numeric=False
    ),
    sr.SymmetryOperation.from_orbitals(
        orbitals=orbitals,
        real_space_operator=sr.RealSpaceOperator(
            rotation_matrix=sp.Matrix([[1, 1, 1], [0, 0, -1], [0, -1, 0]]),
            translation_vector=sp.Matrix([sp.Rational(1, 4)] * 3),
        ),
        rotation_matrix_cartesian=sp.Matrix([[-1, 0, 0], [0, 1, 0], [0, 0,
                                                                     1]]),
        numeric=False
    ),
    sr.SymmetryOperation.from_orbitals(
        orbitals=orbitals,
        real_space_operator=sr.RealSpaceOperator(
            rotation_matrix=sp.Matrix([[1, 0, 0], [-1, -1, -1], [0, 0, 1]]),
            translation_vector=sp.Matrix([0, 0, 0]),
        ),
        rotation_matrix_cartesian=sp.Matrix([[0, 0, -1], [0, 1, 0], [-1, 0,
                                                                     0]]),
        numeric=False
    ),
    sr.get_time_reversal(orbitals=orbitals, numeric=False)
]

print(symmetry_generators)


def print_result(order):
    """prints the basis for a given order of k"""
    print('Order:', order)
    for m in kp.symmetric_hamiltonian(
        *symmetry_generators,
        expr_basis=kp.monomial_basis(order),
        repr_basis=kp.hermitian_basis(16)
    ):
        print(m)
    print()


if __name__ == '__main__':
    for i in range(3):
        print_result(order=i)
