#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    21.01.2017 09:41:11 CET
# File:    matrix_to_vector.py

import sympy

def frobenius_norm(A, B):
    """
    Computes the Frobenius norm <A, B> = Tr(A^\dagger B) for two matrices
    """
    return (A.conjugate().transpose() @ B).trace()
