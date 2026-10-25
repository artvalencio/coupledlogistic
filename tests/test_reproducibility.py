import numpy as np
import pytest

from coupledlogistic import coupledlogistic

XY = [[0, 1], [0, 0]]


def test_same_seed_gives_same_series():
    a = coupledlogistic(300, 4, XY, 0.2, 'diffusive', seed=3)
    b = coupledlogistic(300, 4, XY, 0.2, 'diffusive', seed=3)
    np.testing.assert_array_equal(a, b)


def test_different_seeds_give_different_series():
    a = coupledlogistic(300, 4, XY, 0.2, 'diffusive', seed=3)
    b = coupledlogistic(300, 4, XY, 0.2, 'diffusive', seed=4)
    assert not np.array_equal(a, b)


def test_generator_can_be_passed_as_seed():
    a = coupledlogistic(300, 4, XY, 0.2, 'kaneko', seed=np.random.default_rng(8))
    b = coupledlogistic(300, 4, XY, 0.2, 'kaneko', seed=np.random.default_rng(8))
    np.testing.assert_array_equal(a, b)


def test_raw_output_without_normalization():
    out = coupledlogistic(300, 4, XY, 0.2, 'kaneko', seed=1, normalize=False)
    assert out.shape == (300, 2)
    assert 0 <= out.min() and out.max() <= 1
    assert out.min(axis=0).max() > 0  # not rescaled to [0, 1]


def test_transient_is_discarded():
    full = coupledlogistic(50, 4, XY, 0.2, 'kaneko', seed=2, normalize=False, transient=0)
    late = coupledlogistic(30, 4, XY, 0.2, 'kaneko', seed=2, normalize=False, transient=20)
    np.testing.assert_array_equal(full[20:], late)
