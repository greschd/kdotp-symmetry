#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sympy as sp


def to_matrix(operator, basis, to_vector_fct):
    return sp.Matrix([to_vector_fct(operator(b), basis=basis)
                      for b in basis]).transpose()
