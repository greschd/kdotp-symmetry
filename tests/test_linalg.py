#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from kdotp_symmetry._linalg import zassenhaus, intersection_basis

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

@pytest.mark.parametrize('input_bases', [
    (
        [[0, 1, 0], [0, 1, 1]],
        [[0, 1]],
    ),
    (
        [[0, 1], [0, 1, 1]],
        [[0, 1, 0]],
    ),
    (
        [[0, 1, 0], [0, 1]],
        [[0, 1, 0]],
    ),
])
def test_zassenhaus_inconsistent_dim(input_bases):
    with pytest.raises(ValueError):
            zassenhaus(*input_bases)

@pytest.mark.parametrize('input_bases,output_basis', [
    (
        ([[1, 1, 0], [0, 1, 0]], [[0, 1, 1], [0, 0, 1]], [[0, 2, 1]], []),
        []
    ),
    (
        ([[1, 1, 0], [0, 1, 0]], [[0, 1, 1], [0, 0, 1]]),
        [[0, 1, 0]]
    ),
    (
        ([[1, 1, 0], [0, 1, 0]], [[0, 1, 1], [0, 0, 1]], [(0, 1, 0)]),
        [[0, 1, 0]]
    ),
])
def test_intersection_basis(input_bases, output_basis):
    assert intersection_basis(*input_bases) == output_basis
