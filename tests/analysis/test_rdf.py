import numpy as np
import pytest

from PQAnalysis.analysis.rdf import _calculate_n_bins, _infer_r_max
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


class TestRDF:
    pass
