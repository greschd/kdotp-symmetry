#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    24.01.2017 11:53:31 CET
# File:    test_linalg.py

import pytest
from kdotp_symmetry.linalg import zassenhaus, intersection_basis

@pytest.mark.parametrize('input_bases,output_bases', [
    (
        ([[1, 1, 0], [0, 1, 0]], [[0, 1, 1], [0, 0, 1]]),
        ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 1, 0]])
    ),
    (
        ([[1, 1, 0], [0, 1, 0]], []),
        ([[1, 0, 0], [0, 1, 0]], [])
    ),
    (
        ([[1, 1, 0], [1, 1, 0]], []),
        ([[1, 1, 0]], [])
    ),
    (
        ([[1, 2, 1, 0], [0, -1, 0, 1]], [[1, 1, 0, 1], [0, 0, 1, 0]]),
        ([[1, 0, 0, 2], [0, 1, 0, -1], [0, 0, 1, 0]], [[1, 1, 1, 1]])
    ),
])
def test_zassenhaus(input_bases, output_bases):
    assert zassenhaus(*input_bases) == output_bases
