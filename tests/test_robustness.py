import warnings

import numpy as np
import pytest

from coupledlogistic import coupledlogistic

XY = [[0, 1], [0, 0]]


@pytest.mark.parametrize('seed', range(20))
def test_diverging_orbits_are_restarted(seed):
    # about 1 in 20 of these initial conditions escapes [0, 1] and diverges
    out = coupledlogistic(2000, 4, XY, 0.2, 'diffusive', seed=seed, normalize=False, transient=2000)
    assert np.all(np.isfinite(out))
    assert out.min() >= 0 and out.max() <= 1


def test_impossible_parameters_raise_after_retries():
    with pytest.raises(RuntimeError, match='diverg'):
        coupledlogistic(100, 4.5, XY, 0.2, 'diffusive', seed=0, max_retries=3)


def test_retries_are_silent_unless_verbose(capsys):
    with pytest.raises(RuntimeError):
        coupledlogistic(100, 4.5, XY, 0.2, 'diffusive', seed=0, max_retries=2)
    assert capsys.readouterr().out == ''
    with pytest.raises(RuntimeError):
        coupledlogistic(100, 4.5, XY, 0.2, 'diffusive', seed=0, max_retries=2, verbose=True)
    assert 'recalculating' in capsys.readouterr().out


def test_constant_series_normalize_without_nan():
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        out = coupledlogistic(1, 4, XY, 0.2, 'kaneko', seed=0)
    np.testing.assert_array_equal(out, 0)


def test_verbose_works_for_short_series(capsys):
    coupledlogistic(5, 4, XY, 0.2, 'kaneko', seed=0, transient=0, verbose=True)
    assert 'Generating' in capsys.readouterr().out
