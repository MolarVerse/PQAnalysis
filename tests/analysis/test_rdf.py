import numpy as np
import pytest

from PQAnalysis.analysis.rdf import _calculate_n_bins, _infer_r_max, _integration, _norm
from PQAnalysis.analysis import RDFError
from PQAnalysis.traj import Trajectory, Frame
from PQAnalysis.core import AtomicSystem, Cell


def test__calculate_n_bins():
    r_min = 1.0
    r_max = 101.5
    delta_r = 1.0

    n_bins, r_max = _calculate_n_bins(r_min, r_max, delta_r)

    assert n_bins == 100
    assert np.isclose(r_max, 101.0)


def test__infer_r_max():

    system1 = AtomicSystem(cell=Cell(10, 10, 10, 90, 90, 90))
    system2 = AtomicSystem(cell=Cell(16, 13, 12, 90, 90, 90))

    traj = Trajectory([Frame(system1), Frame(system2)])

    r_max = _infer_r_max(traj)

    assert np.isclose(r_max, 5.0)

    system3 = AtomicSystem()
    traj.append(Frame(system3))

    print(traj.check_vacuum())

    with pytest.raises(RDFError) as exception:
        r_max = _infer_r_max(traj)
    assert str(exception.value) == "To infer r_max of the RDF analysis, the trajectory cannot be a vacuum trajectory. Please specify r_max manually or use the combination n_bins and delta_r."


def test__integration():
    bins = np.array([1, 2, 3, 4, 5])
    len_reference_indices = 3
    len_frames = 10

    integration = _integration(bins, len_reference_indices, len_frames)

    n_total = len_reference_indices * len_frames
    assert np.allclose(integration, np.array(
        [1 / n_total, 3 / n_total, 6 / n_total, 10 / n_total, 15 / n_total]))


def test__norm():
    n_bins = 5
    n_frames = 10
    n_reference_indices = 3
    delta_r = 1.0
    target_density = 2.0

    norm = _norm(n_bins, delta_r, target_density,
                 n_reference_indices, n_frames)

    help_1 = np.arange(0, n_bins)
    help_2 = np.arange(1, n_bins + 1)
    norm_ref = (help_2**3 - help_1**3)*delta_r**3 * 4 / 3 * np.pi

    assert np.allclose(norm, norm_ref * target_density *
                       n_reference_indices * n_frames)


class TestRDF:
    pass
