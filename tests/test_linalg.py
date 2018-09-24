"""
Tests for the Zassenhaus algorithm and basis intersection.
"""

import pytest
from kdotp_symmetry._linalg import zassenhaus, intersection_basis


@pytest.mark.parametrize(
    'input_bases,output_bases', [
        (([[1, 1, 0], [0, 1, 0]], [[0, 1, 1], [0, 0, 1]]),
         ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 1, 0]])),
        (([[1, 1, 0], [0, 1, 0]], []), ([[1, 0, 0], [0, 1, 0]], [])),
        (([[1, 1, 0], [1, 1, 0]], []), ([[1, 1, 0]], [])),
        (([[1, 2, 1, 0], [0, -1, 0, 1]], [[1, 1, 0, 1], [0, 0, 1, 0]]),
         ([[1, 0, 0, 2], [0, 1, 0, -1], [0, 0, 1, 0]], [[1, 1, 1, 1]])),
    ]
)
def test_zassenhaus(input_bases, output_bases):
    """
    Test the Zassenhaus algorithm.
    """
    assert zassenhaus(*input_bases) == output_bases


@pytest.mark.parametrize(
    'input_bases', [
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
    ]
)
def test_zassenhaus_inconsistent_dim(input_bases):
    """
    Test that the Zassenhaus algorithm raises an error when the input has inconsistent dimension.
    """
    with pytest.raises(ValueError):
        zassenhaus(*input_bases)


@pytest.mark.parametrize(
    'input_bases,output_basis', [
        (([[1, 1, 0], [0, 1, 0]], [[0, 1, 1], [0, 0, 1]], [[0, 2, 1]], []),
         []),
        (([[1, 1, 0], [0, 1, 0]], [[0, 1, 1], [0, 0, 1]]), [[0, 1, 0]]),
        (([[1, 1, 0], [0, 1, 0]], [[0, 1, 1], [0, 0, 1]], [(0, 1, 0)]),
         [[0, 1, 0]]),
    ]
)
def test_intersection_basis(input_bases, output_basis):
    """
    Test the basis intersection algorithm.
    """
    assert intersection_basis(*input_bases) == output_basis
