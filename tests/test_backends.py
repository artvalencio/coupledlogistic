import importlib.util
import time

import numpy as np
import pytest

from coupledlogistic import coupledlogistic

HAS_NUMBA = importlib.util.find_spec('numba') is not None
WHEATSTONE = [[0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 0, 1], [0, 0, 0, 0]]
needs_numba = pytest.mark.skipif(not HAS_NUMBA, reason='numba not installed')


@needs_numba
@pytest.mark.parametrize('couplingtype', ['diffusive', 'kaneko'])
def test_backends_agree_before_chaos_amplifies_rounding(couplingtype):
    # rounding differences grow by ~2x per step (r=4), so only the first steps can be compared
    kw = dict(seed=7, normalize=False, transient=0)
    a = coupledlogistic(10, 4, WHEATSTONE, 0.15, couplingtype, backend='numpy', **kw)
    b = coupledlogistic(10, 4, WHEATSTONE, 0.15, couplingtype, backend='numba', **kw)
    np.testing.assert_allclose(a, b, rtol=0, atol=1e-12)


@pytest.mark.parametrize('backend', ['numpy', pytest.param('numba', marks=needs_numba)])
def test_seed_reproduces_within_a_backend(backend):
    a = coupledlogistic(2000, 4, WHEATSTONE, 0.15, 'kaneko', seed=3, backend=backend)
    b = coupledlogistic(2000, 4, WHEATSTONE, 0.15, 'kaneko', seed=3, backend=backend)
    np.testing.assert_array_equal(a, b)


def test_unknown_backend_is_rejected():
    with pytest.raises(ValueError, match='backend'):
        coupledlogistic(10, 4, WHEATSTONE, 0.15, 'kaneko', backend='fortran')


@pytest.mark.parametrize('backend', ['numpy', pytest.param('numba', marks=needs_numba)])
def test_divergence_is_detected_early(backend):
    # each failed attempt must stop at the divergence, not simulate the full length
    start = time.perf_counter()
    with pytest.raises(RuntimeError, match='diverg'):
        coupledlogistic(10**6, 4.5, [[0, 1], [0, 0]], 0.2, 'diffusive', seed=0, backend=backend)
    assert time.perf_counter() - start < 1.0


@pytest.mark.parametrize('backend', ['numpy', pytest.param('numba', marks=needs_numba)])
def test_float32_output(backend):
    out64 = coupledlogistic(500, 4, WHEATSTONE, 0.15, 'kaneko', seed=1, backend=backend)
    out32 = coupledlogistic(500, 4, WHEATSTONE, 0.15, 'kaneko', seed=1, backend=backend, dtype=np.float32)
    assert out32.dtype == np.float32
    np.testing.assert_allclose(out32, out64, atol=1e-6)


@pytest.mark.parametrize('backend', ['numpy', pytest.param('numba', marks=needs_numba)])
def test_realizations(backend):
    from test_dynamics import assert_follows_equations
    out = coupledlogistic(300, 4, WHEATSTONE, 0.15, 'diffusive', seed=2, normalize=False, transient=100,
                          n_realizations=5, backend=backend)
    assert out.shape == (5, 300, 4)
    for realization in out:
        assert_follows_equations(realization, 4, WHEATSTONE, 0.15, 'diffusive')
    assert not np.allclose(out[0], out[1])
    again = coupledlogistic(300, 4, WHEATSTONE, 0.15, 'diffusive', seed=2, normalize=False, transient=100,
                            n_realizations=5, backend=backend)
    np.testing.assert_array_equal(out, again)


def test_normalized_realizations_are_scaled_separately():
    out = coupledlogistic(500, 4, WHEATSTONE, 0.15, 'kaneko', seed=2, n_realizations=3)
    np.testing.assert_allclose(out.min(axis=1), 0)
    np.testing.assert_allclose(out.max(axis=1), 1)
