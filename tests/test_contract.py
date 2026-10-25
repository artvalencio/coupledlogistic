import numpy as np
import pytest

from coupledlogistic import coupledlogistic

XY = [[0, 1], [0, 0]]  # node 0 drives node 1


@pytest.mark.parametrize('couplingtype', ['diffusive', 'kaneko'])
@pytest.mark.parametrize('n', [1, 10, 500])
def test_output_has_requested_shape(couplingtype, n):
    assert coupledlogistic(n, 4, XY, 0.2, couplingtype).shape == (n, 2)


def test_output_is_normalized_to_unit_interval():
    out = coupledlogistic(1000, 4, XY, 0.2, 'diffusive')
    np.testing.assert_allclose(out.min(axis=0), 0)
    np.testing.assert_allclose(out.max(axis=0), 1)


def test_unknown_coupling_type_is_rejected():
    with pytest.raises(ValueError, match='couplingtype'):
        coupledlogistic(50, 4, XY, 0.2, 'magnetic')


def test_adjacency_matrix_must_be_square():
    with pytest.raises(ValueError, match='square'):
        coupledlogistic(50, 4, [[0, 1, 0], [0, 0, 1]], 0.2, 'diffusive')


def test_r_vector_must_match_number_of_nodes():
    with pytest.raises(ValueError, match='r'):
        coupledlogistic(50, [4, 3.9, 3.8], XY, 0.2, 'diffusive')


def test_single_node_network():
    assert coupledlogistic(100, 4, [[0]], 0.2, 'diffusive').shape == (100, 1)
