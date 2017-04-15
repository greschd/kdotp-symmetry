#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    21.01.2017 20:45:51 CET
# File:    to_matrix.py

import sympy as sp

def to_matrix(operator, basis, to_vector_fct):
    return  sp.Matrix([
        to_vector_fct(operator(b), basis=basis)
        for b in basis
    ]).transpose()
