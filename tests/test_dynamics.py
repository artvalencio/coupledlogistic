"""The simulated orbits must follow the equations in the docstring exactly.

diffusive: x_i' = (1-s) f_i(x_i) + s/k_i sum_j A_ji (x_j - x_i)
kaneko:    x_i' = (1-s) f_i(x_i) + s/k_i sum_j A_ji f_j(x_j)
with f_i(x) = r_i x (1-x), k_i = sum_j A_ji, and x_i' = f_i(x_i) when k_i = 0.
"""
import numpy as np
import pytest

from coupledlogistic import coupledlogistic

NETWORKS = {
    'pair': [[0, 1], [0, 0]],
    'serial': [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0]],
    'parallel': [[0, 1, 1, 0], [0, 0, 0, 1], [0, 0, 0, 1], [0, 0, 0, 0]],
    'wheatstone': [[0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 0, 1], [0, 0, 0, 0]],
}


def reference_step(x, r, A, sigma, couplingtype):
    A = np.asarray(A, dtype=float)
    f = r * x * (1 - x)
    new = np.empty_like(x)
    for i in range(len(x)):
        k = sum(A[j, i] for j in range(len(x)))
        if k == 0:
            new[i] = f[i]
            continue
        if couplingtype == 'diffusive':
            coupling = sum(A[j, i] * (x[j] - x[i]) for j in range(len(x)))
        else:
            coupling = sum(A[j, i] * f[j] for j in range(len(x)))
        new[i] = (1 - sigma) * f[i] + sigma * coupling / k
    return new


def assert_follows_equations(out, r, A, sigma, couplingtype):
    r = np.broadcast_to(np.asarray(r, dtype=float), out.shape[1])
    expected = np.array([reference_step(row, r, A, sigma, couplingtype) for row in out[:-1]])
    np.testing.assert_allclose(out[1:], expected, rtol=1e-12, atol=1e-15)


BACKENDS = ['numpy', pytest.param('numba', marks=pytest.mark.skipif(
    __import__('importlib').util.find_spec('numba') is None, reason='numba not installed'))]


@pytest.mark.parametrize('backend', BACKENDS)
@pytest.mark.parametrize('couplingtype', ['diffusive', 'kaneko'])
@pytest.mark.parametrize('network', NETWORKS)
def test_orbit_follows_equations(network, couplingtype, backend):
    A = NETWORKS[network]
    out = coupledlogistic(200, 4, A, 0.15, couplingtype, seed=11, normalize=False, transient=0, backend=backend)
    assert_follows_equations(out, 4, A, 0.15, couplingtype)


@pytest.mark.parametrize('backend', BACKENDS)
@pytest.mark.parametrize('couplingtype', ['diffusive', 'kaneko'])
def test_orbit_after_transient_follows_equations(couplingtype, backend):
    A = NETWORKS['wheatstone']
    out = coupledlogistic(300, [3.9, 3.7, 4.0, 3.8], A, 0.15, couplingtype, seed=4, normalize=False, transient=5000,
                          backend=backend)
    assert_follows_equations(out, [3.9, 3.7, 4.0, 3.8], A, 0.15, couplingtype)


@pytest.mark.parametrize('couplingtype', ['diffusive', 'kaneko'])
def test_heterogeneous_r(couplingtype):
    A, r = NETWORKS['parallel'], [3.9, 3.7, 4.0, 3.8]
    out = coupledlogistic(200, r, A, 0.15, couplingtype, seed=5, normalize=False, transient=0)
    assert_follows_equations(out, r, A, 0.15, couplingtype)


@pytest.mark.parametrize('couplingtype', ['diffusive', 'kaneko'])
def test_weighted_adjacency(couplingtype):
    A = [[0, 0.5, 2.0], [0, 0, 1.0], [0, 0, 0]]
    out = coupledlogistic(200, 4, A, 0.1, couplingtype, seed=2, normalize=False, transient=0)
    assert_follows_equations(out, 4, A, 0.1, couplingtype)


def test_without_coupling_every_node_is_a_logistic_map():
    out = coupledlogistic(100, [4, 3.9], NETWORKS['pair'], 0.0, 'kaneko', seed=1, normalize=False, transient=0)
    np.testing.assert_allclose(out[1:], np.array([4, 3.9]) * out[:-1] * (1 - out[:-1]))


@pytest.mark.parametrize('backend', BACKENDS)
def test_random_weighted_network(backend):
    rng = np.random.default_rng(0)
    A = (rng.random((12, 12)) < 0.2) * rng.random((12, 12))
    np.fill_diagonal(A, 0)
    for couplingtype in ('diffusive', 'kaneko'):
        out = coupledlogistic(150, 4, A, 0.1, couplingtype, seed=3, normalize=False, transient=0, backend=backend)
        assert_follows_equations(out, 4, A, 0.1, couplingtype)
