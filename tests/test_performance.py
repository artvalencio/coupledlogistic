import time

import numpy as np

import pytest

from coupledlogistic import coupledlogistic

WHEATSTONE = [[0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 0, 1], [0, 0, 0, 0]]


@pytest.mark.slow
@pytest.mark.parametrize('couplingtype', ['diffusive', 'kaneko'])
def test_long_series_are_generated_quickly(couplingtype):
    start = time.perf_counter()
    coupledlogistic(100_000, 4, WHEATSTONE, 0.15, couplingtype, seed=1)  # + 10000 transient steps
    assert time.perf_counter() - start < 2.0


@pytest.mark.slow
def test_numba_backend_is_fast():
    pytest.importorskip('numba')
    coupledlogistic(10, 4, WHEATSTONE, 0.15, 'kaneko', seed=1, backend='numba')  # compile
    start = time.perf_counter()
    coupledlogistic(10**6, 4, WHEATSTONE, 0.15, 'diffusive', seed=1, backend='numba')
    assert time.perf_counter() - start < 0.2


@pytest.mark.slow
def test_many_realizations_numpy_backend():
    start = time.perf_counter()
    coupledlogistic(10**4, 4, WHEATSTONE, 0.15, 'kaneko', seed=1, transient=1000, n_realizations=100, backend='numpy')
    assert time.perf_counter() - start < 1.0


@pytest.mark.slow
@pytest.mark.parametrize('dtype,limit', [(np.float64, 8), (np.float32, 4.5)])
def test_memory(dtype, limit):
    import tracemalloc
    tracemalloc.start()
    coupledlogistic(2 * 10**5, 4, WHEATSTONE, 0.15, 'kaneko', seed=1, backend='numpy', dtype=dtype, transient=1000)
    peak = tracemalloc.get_traced_memory()[1] / 1e6
    tracemalloc.stop()
    assert peak < limit  # output alone: 6.4 MB (float64) or 3.2 MB (float32)


@pytest.mark.slow
def test_numba_large_sparse_network_is_fast():
    pytest.importorskip('numba')
    chain = np.eye(50, k=1)
    coupledlogistic(10, 4, chain, 0.15, 'kaneko', seed=1, backend='numba')  # compile
    start = time.perf_counter()
    coupledlogistic(10**5, 4, chain, 0.15, 'kaneko', seed=1, backend='numba')
    assert time.perf_counter() - start < 0.06
