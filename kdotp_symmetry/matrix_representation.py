#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    21.01.2017 20:45:51 CET
# File:    matrix_representation.py

import sympy as sp

def matrix_representation(operator, basis, to_vector_fct):
    print([operator(b) for b in basis])
    return sp.Matrix([
        to_vector_fct(operator(b), basis=basis)
        for b in basis
    ]).transpose()
