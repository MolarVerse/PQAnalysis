import numpy as np
import pytest

from PQAnalysis.analysis.rdf.exceptions import RDFError
from PQAnalysis.analysis import RDF
from PQAnalysis.traj import Trajectory
from PQAnalysis.core import Cell
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.topology import SelectionCompatible
from PQAnalysis.types import PositiveReal, PositiveInt
from PQAnalysis.exceptions import PQTypeError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception

# pylint: disable=protected-access



def test__calculate_n_bins():
    r_min = 1.0
    r_max = 101.5
    delta_r = 1.0

    n_bins, r_max = RDF._calculate_n_bins(r_min, r_max, delta_r)

    assert n_bins == 100
    assert np.isclose(r_max, 101.0)



def test__infer_r_max(caplog):

    system1 = AtomicSystem(cell=Cell(10, 10, 10, 90, 90, 90))
    system2 = AtomicSystem(cell=Cell(16, 13, 12, 90, 90, 90))

    traj = Trajectory([system1, system2])

    r_max = RDF._infer_r_max(traj.cells)

    assert np.isclose(r_max, 5.0)

    system3 = AtomicSystem()
    traj.append(system3)

    assert_logging_with_exception(
        caplog=caplog,
        logging_name=RDF.__qualname__,
        logging_level="ERROR",
        message_to_test=(
        "To infer r_max of the RDF analysis, the "
        "trajectory cannot be a vacuum trajectory. "
        "Please specify r_max manually or use "
        "the combination n_bins and delta_r."
        ),
        exception=RDFError,
        function=RDF._infer_r_max,
        cells=traj.cells,
    )



def test__check_r_max(caplog):
    r_max = 5.0
    traj = Trajectory()

    assert np.isclose(RDF._check_r_max(r_max, traj.cells), r_max)

    system1 = AtomicSystem(cell=Cell(10, 10, 10, 90, 90, 90))
    system2 = AtomicSystem(cell=Cell(16, 13, 12, 90, 90, 90))

    traj = Trajectory([system1, system2])

    assert np.isclose(RDF._check_r_max(r_max, traj.cells), r_max)

    r_max = 10.0

    r_max = assert_logging_with_exception(
        caplog=caplog,
        logging_name=RDF.__qualname__,
        logging_level="WARNING",
        message_to_test=(
        f"The calculated r_max {r_max} "
        "is larger than the maximum allowed radius "
        "according to the box vectors of the trajectory 5.0. "
        "r_max will be set to the maximum allowed radius."
        ),
        exception=None,
        function=RDF._check_r_max,
        r_max=r_max,
        cells=traj.cells,
    )

    assert np.isclose(r_max, 5.0)



def test__calculate_r_max(caplog):
    n_bins = 50
    delta_r = 0.1
    r_min = 0.0
    traj = Trajectory()

    r_max = RDF._calculate_r_max(n_bins, delta_r, r_min, traj.cells)

    assert np.isclose(r_max, 5.0)

    r_min = 3.0
    r_max = RDF._calculate_r_max(n_bins, delta_r, r_min, traj.cells)

    assert np.isclose(r_max, 8.0)



def test__setup_bin_middle_points():
    n_bins = 5
    r_min = 3.0
    r_max = 8.0
    delta_r = 1.0

    bin_middle_points = RDF._setup_bin_middle_points(
        n_bins,
        r_min,
        r_max,
        delta_r
    )

    assert np.allclose(bin_middle_points, np.array([3.5, 4.5, 5.5, 6.5, 7.5]))



def test__integration():
    bins = np.array([1, 2, 3, 4, 5])
    len_reference_indices = 3
    len_frames = 10

    integration = RDF._integration(bins, len_reference_indices, len_frames)

    n_total = len_reference_indices * len_frames
    assert np.allclose(
        integration,
        np.array(
        [1 / n_total,
        3 / n_total,
        6 / n_total,
        10 / n_total,
        15 / n_total]
        )
    )



def test__norm():
    n_bins = 5
    n_frames = 10
    n_reference_indices = 3
    delta_r = 1.0
    target_density = 2.0

    norm = RDF._norm(
        n_bins,
        delta_r,
        target_density,
        n_reference_indices,
        n_frames
    )

    help_1 = np.arange(0, n_bins)
    help_2 = np.arange(1, n_bins + 1)
    norm_ref = (help_2**3 - help_1**3) * delta_r**3 * 4 / 3 * np.pi

    assert np.allclose(
        norm,
        norm_ref * target_density * n_reference_indices * n_frames
    )



def test__add_to_bins():
    n_bins = 5
    r_min = 3.0
    delta_r = 1.0

    distances = np.array([1.5, 2.5, 3.5, 3.6, 3.7, 4.5, 4.6, 5.5, 6.5, 8.5])

    assert np.allclose(
        RDF._add_to_bins(distances,
        r_min,
        delta_r,
        n_bins),
        np.array([3,
        2,
        1,
        1,
        0])
    )



class TestRDF:

    def test__init__type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "traj",
            1,
            Trajectory | TrajectoryReader
            ),
            exception=PQTypeError,
            function=RDF,
            traj=1,
            reference_species=["h"],
            target_species=["h"],
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "reference_species",
            Trajectory(),
            SelectionCompatible
            ),
            exception=PQTypeError,
            function=RDF,
            traj=Trajectory(),
            reference_species=Trajectory(),
            target_species=["h"],
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "target_species",
            Trajectory(),
            SelectionCompatible
            ),
            exception=PQTypeError,
            function=RDF,
            traj=Trajectory(),
            reference_species=["h"],
            target_species=Trajectory(),
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "use_full_atom_info",
            1,
            bool
            ),
            exception=PQTypeError,
            function=RDF,
            traj=Trajectory(),
            reference_species=["h"],
            target_species=["h"],
            use_full_atom_info=1,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "no_intra_molecular",
            1,
            bool
            ),
            exception=PQTypeError,
            function=RDF,
            traj=Trajectory(),
            reference_species=["h"],
            target_species=["h"],
            no_intra_molecular=1,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "r_max",
            -1,
            PositiveReal | None
            ),
            exception=PQTypeError,
            function=RDF,
            traj=Trajectory(),
            reference_species=["h"],
            target_species=["h"],
            r_max=-1,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "r_min",
            -1,
            PositiveReal | None
            ),
            exception=PQTypeError,
            function=RDF,
            traj=Trajectory(),
            reference_species=["h"],
            target_species=["h"],
            r_min=-1,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "delta_r",
            -1,
            PositiveReal | None
            ),
            exception=PQTypeError,
            function=RDF,
            traj=Trajectory(),
            reference_species=["h"],
            target_species=["h"],
            delta_r=-1,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "n_bins",
            -1,
            PositiveInt | None
            ),
            exception=PQTypeError,
            function=RDF,
            traj=Trajectory(),
            reference_species=["h"],
            target_species=["h"],
            n_bins=-1,
        )

    def test__init__(self, caplog):

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test="Trajectory cannot be of length 0.",
            exception=RDFError,
            function=RDF,
            traj=Trajectory(),
            reference_species=["h"],
            target_species=["h"],
            r_max=8.0,
            r_min=3.0,
        )

        system1 = AtomicSystem(cell=Cell(10, 10, 10, 90, 90, 90))
        system2 = AtomicSystem(cell=Cell())
        traj = Trajectory([system1, system2])

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
            "The provided trajectory is not fully periodic "
            "or in vacuum, meaning that some frames are in "
            "vacuum and others are periodic. This is not "
            "supported by the RDF analysis."
            ),
            exception=RDFError,
            function=RDF,
            traj=traj,
            reference_species=["h"],
            target_species=["h"],
            r_max=8.0,
            r_min=3.0,
        )

        system1 = AtomicSystem(cell=Cell(10, 10, 10, 90, 90, 90))
        system2 = AtomicSystem(cell=Cell(16, 13, 12, 90, 90, 90))

        traj = Trajectory([system1, system2])

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test="Either n_bins or delta_r must be specified.",
            exception=RDFError,
            function=RDF,
            traj=traj,
            reference_species=["h"],
            target_species=["h"],
            r_max=8.0,
            r_min=3.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
            "It is not possible to specify all of n_bins, "
            "delta_r and r_max in the same RDF analysis as "
            "this would lead to ambiguous results."
            ),
            exception=RDFError,
            function=RDF,
            traj=traj,
            reference_species=["h"],
            target_species=["h"],
            r_max=8.0,
            r_min=3.0,
            delta_r=0.1,
            n_bins=5,
        )

        # initialize rdf only with n_bins and delta_r

        n_bins = 5
        delta_r = 1.0

        rdf = RDF(traj, ["h"], ["h"], delta_r=delta_r, n_bins=n_bins)

        assert np.isclose(rdf.r_max, 5.0)
        assert np.isclose(rdf.r_min, 0.0)
        assert len(rdf.bins) == 5
        assert np.allclose(
            rdf.bin_middle_points,
            np.array([0.5,
            1.5,
            2.5,
            3.5,
            4.5])
        )
        assert rdf.n_bins == 5
        assert np.isclose(rdf.delta_r, 1.0)

        # r_max has to be taken from trajectory

        n_bins = 10

        rdf = assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="WARNING",
            message_to_test=(
            "The calculated r_max 10.0 is larger than the maximum allowed "
            "radius according to the box vectors of the trajectory 5.0. "
            "r_max will be set to the maximum allowed radius."
            ),
            exception=None,
            function=RDF,
            traj=traj,
            reference_species=["h"],
            target_species=["h"],
            delta_r=delta_r,
            n_bins=n_bins,
        )

        assert np.isclose(rdf.r_max, 5.0)

        system1 = AtomicSystem(cell=Cell())
        system2 = AtomicSystem(cell=Cell())

        traj = Trajectory([system1, system2])

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
            "To infer r_max of the RDF analysis, the trajectory cannot "
            "be a vacuum trajectory. Please specify r_max manually or "
            "use the combination n_bins and delta_r."
            ),
            exception=RDFError,
            function=RDF,
            traj=traj,
            reference_species=["h"],
            target_species=["h"],
            delta_r=delta_r,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
            "To infer r_max of the RDF analysis, the trajectory cannot be "
            "a vacuum trajectory. Please specify r_max manually or use the "
            "combination n_bins and delta_r."
            ),
            exception=RDFError,
            function=RDF,
            traj=traj,
            reference_species=["h"],
            target_species=["h"],
            n_bins=n_bins,
        )

        r_max = 5.0

        rdf = RDF(traj, ["h"], ["h"], delta_r=delta_r, r_max=r_max)

        assert np.isclose(rdf.r_max, 5.0)
        assert np.isclose(rdf.r_min, 0.0)
        assert len(rdf.bins) == 5
        assert np.allclose(
            rdf.bin_middle_points,
            np.array([0.5,
            1.5,
            2.5,
            3.5,
            4.5])
        )
        assert rdf.n_bins == 5
        assert np.isclose(rdf.delta_r, 1.0)

        n_bins = 5

        rdf = RDF(traj, ["h"], ["h"], n_bins=n_bins, r_max=r_max)

        assert np.isclose(rdf.r_max, 5.0)
        assert np.isclose(rdf.r_min, 0.0)
        assert len(rdf.bins) == 5
        assert np.allclose(
            rdf.bin_middle_points,
            np.array([0.5,
            1.5,
            2.5,
            3.5,
            4.5])
        )
        assert rdf.n_bins == 5
        assert np.isclose(rdf.delta_r, 1.0)
