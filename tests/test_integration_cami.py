import pytest

from coupledlogistic import coupledlogistic

cami = pytest.importorskip('cami')


@pytest.mark.parametrize('couplingtype', ['diffusive', 'kaneko'])
def test_information_flows_along_the_link(couplingtype):
    out = coupledlogistic(20000, 4, [[0, 1], [0, 0]], 0.2, couplingtype, seed=0)
    assert cami.directionality(out[:, 0], out[:, 1], symbolic_length=2) > 0.05
