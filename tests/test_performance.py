import time

import pytest

from coupledlogistic import coupledlogistic

WHEATSTONE = [[0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 0, 1], [0, 0, 0, 0]]


@pytest.mark.slow
@pytest.mark.parametrize('couplingtype', ['diffusive', 'kaneko'])
def test_long_series_are_generated_quickly(couplingtype):
    start = time.perf_counter()
    coupledlogistic(100_000, 4, WHEATSTONE, 0.15, couplingtype, seed=1)  # + 10000 transient steps
    assert time.perf_counter() - start < 2.0
