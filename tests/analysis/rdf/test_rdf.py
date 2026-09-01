import sys

import numpy as np
import pytest

from PQAnalysis import config
from PQAnalysis.analysis.rdf.exceptions import RDFError
from PQAnalysis.analysis import RDF
from PQAnalysis.traj import Trajectory
from PQAnalysis.core import Atom, Cell, Element, Residue
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.topology import SelectionCompatible, Topology
from PQAnalysis.types import PositiveReal, PositiveInt
from PQAnalysis.exceptions import PQTypeError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception

# pylint: disable=protected-access



def _make_no_intra_trajectory():
    system1 = AtomicSystem(
        atoms=[Atom("H"), Atom("H"), Atom("C")],
        pos=np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]]),
        cell=Cell(10, 10, 10, 90, 90, 90)
    )
    system2 = AtomicSystem(
        atoms=[Atom("H"), Atom("H"), Atom("C")],
        pos=np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]]),
        cell=Cell(10, 10, 10, 90, 90, 90)
    )

    return Trajectory([system1, system2])


def _make_residue_trajectory(n_water: int, n_sodium: int):
    """
    Builds a trajectory whose topology has real residues:
    ``n_water`` water molecules (3 atoms each) followed by
    ``n_sodium`` single-atom sodium residues.
    """
    water = Residue(
        name="WAT",
        residue_id=1,
        total_charge=0.0,
        elements=[Element("O"), Element("H"), Element("H")],
        atom_types=np.array([0, 1, 1]),
        partial_charges=np.array([-0.8, 0.4, 0.4]),
    )

    atoms = [Atom("O"), Atom("H"), Atom("H")] * n_water
    atoms += [Atom("Na")] * n_sodium
    residue_ids = np.array([1] * (3 * n_water) + [0] * n_sodium)

    topology = Topology(
        atoms=atoms,
        residue_ids=residue_ids,
        reference_residues=[water],
    )

    positions = []
    for i in range(n_water):
        origin = 4.3 * i
        positions.append([origin, 0.0, 0.0])
        positions.append([origin + 1.0, 0.0, 0.0])
        positions.append([origin, 1.0, 0.0])
    for i in range(n_sodium):
        positions.append([float(i), 12.0, 12.0])

    system = AtomicSystem(
        pos=np.array(positions),
        cell=Cell(30, 30, 30, 90, 90, 90),
        topology=topology,
    )

    return Trajectory([system])


def _make_partial_rdf_reference_trajectory():
    symbols = ["H", "H", "H", "O", "O", "O", "O"]
    box_length = 12.0
    cell = Cell(box_length, box_length, box_length, 90, 90, 90)
    positions_by_frame = [
        np.array(
            [
                [0.4, 0.7, 1.1],
                [5.3, 5.5, 5.1],
                [10.8, 1.2, 4.4],
                [1.2, 0.9, 1.4],
                [6.7, 5.1, 5.8],
                [11.5, 11.7, 4.9],
                [3.1, 7.4, 10.6],
            ]
        ),
        np.array(
            [
                [0.6, 0.9, 1.0],
                [5.0, 5.8, 5.2],
                [10.6, 1.4, 4.2],
                [1.5, 1.1, 1.8],
                [6.4, 5.4, 5.7],
                [11.3, 11.6, 5.2],
                [3.4, 7.0, 10.1],
            ]
        ),
        np.array(
            [
                [0.2, 0.5, 1.4],
                [5.6, 5.2, 4.8],
                [10.9, 1.0, 4.6],
                [1.0, 1.3, 1.7],
                [6.9, 4.8, 5.4],
                [11.6, 11.4, 4.6],
                [3.3, 7.7, 10.9],
            ]
        ),
    ]
    trajectory = Trajectory([
        AtomicSystem(
            atoms=[Atom(symbol) for symbol in symbols],
            pos=positions,
            cell=cell,
        )
        for positions in positions_by_frame
    ])

    return symbols, box_length, positions_by_frame, trajectory



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



@pytest.mark.parametrize("r_min", [0.0, 3.0])
def test__norm(r_min):
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
        n_frames,
        r_min,
    )

    inner_radii = r_min + np.arange(n_bins) * delta_r
    outer_radii = inner_radii + delta_r
    norm_ref = (outer_radii**3 - inner_radii**3) * 4 / 3 * np.pi

    assert np.allclose(
        norm,
        norm_ref * target_density * n_reference_indices * n_frames
    )



def test_run_normalizes_shells_from_nonzero_r_min():
    cell = Cell(10.0, 10.0, 10.0)
    trajectory = Trajectory([
        AtomicSystem(
            atoms=[Atom("H"), Atom("O")],
            pos=np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]),
            cell=cell,
        )
    ])

    centers, rdf_values, coordination, shell_population, residual = RDF(
        trajectory,
        "H",
        "O",
        n_bins=1,
        delta_r=1.0,
        r_min=1.0,
    ).run()

    shell_volume = 4.0 / 3.0 * np.pi * (2.0**3 - 1.0**3)
    expected_count = shell_volume / cell.volume

    assert np.allclose(centers, [1.5])
    assert np.allclose(rdf_values, [1.0 / expected_count])
    assert np.allclose(coordination, [1.0])
    assert np.allclose(shell_population, [cell.volume])
    assert np.allclose(residual, [1.0 - expected_count])



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



def test_progress_bar_binds_config_at_call_time(monkeypatch):
    # config.with_progress_bar is set by the CLI after the module
    # import, so it must be read at call time, not bound by value
    # at import time
    captured = {}

    def fake_tqdm(iterable, **kwargs):
        captured.update(kwargs)
        return iterable

    rdf_module = sys.modules[RDF.__module__]
    monkeypatch.setattr(rdf_module, "tqdm", fake_tqdm)

    monkeypatch.setattr(config, "with_progress_bar", False)
    RDF(_make_no_intra_trajectory(), ["H"], ["C"], delta_r=0.5, n_bins=5).run()
    assert captured["disable"] is True

    captured.clear()

    monkeypatch.setattr(config, "with_progress_bar", True)
    RDF(_make_no_intra_trajectory(), ["H"], ["C"], delta_r=0.5, n_bins=5).run()
    assert captured["disable"] is False



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
            bool | None
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
            bool | None
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

        system1 = AtomicSystem(
            atoms=[Atom("h")],
            pos=np.array([[0, 0, 0]]),
            cell=Cell(10, 10, 10, 90, 90, 90)
        )
        system2 = AtomicSystem(
            atoms=[Atom("h")],
            pos=np.array([[0, 0, 0]]),
            cell=Cell(16, 13, 12, 90, 90, 90)
        )

        traj = Trajectory([system1, system2])

        rdf = RDF(
            traj=traj,
            reference_species=["h"],
            target_species=["h"],
            delta_r=0.1,
            use_full_atom_info=None,
            no_intra_molecular=None,
            r_min=None,
        )

        assert rdf.use_full_atom_info == rdf._use_full_atom_default
        assert rdf.no_intra_molecular == rdf._no_intra_molecular_default
        assert rdf.r_min == rdf._r_min_default

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

        system1 = AtomicSystem(
            atoms=[Atom("h")], pos=np.array([[0, 0, 0]]), cell=Cell()
        )
        system2 = AtomicSystem(
            atoms=[Atom("h")], pos=np.array([[0, 0, 0]]), cell=Cell()
        )

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

    def test__init__empty_selection(self, caplog):
        system = AtomicSystem(
            atoms=[Atom("O"), Atom("H")],
            pos=np.array([[0, 0, 0], [1, 0, 0]]),
            cell=Cell(10, 10, 10, 90, 90, 90)
        )
        traj = Trajectory([system])

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
            "The reference selection does not select any atoms."
            ),
            exception=RDFError,
            function=RDF,
            traj=traj,
            reference_species=["Na"],
            target_species=["H"],
            delta_r=0.5,
            r_max=4.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
            "The target selection does not select any atoms."
            ),
            exception=RDFError,
            function=RDF,
            traj=traj,
            reference_species=["O"],
            target_species=["Na"],
            delta_r=0.5,
            r_max=4.0,
        )

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test__init__empty_selection_legacy_path(
        self, caplog, test_with_data_dir
    ):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
            "The reference selection does not select any atoms."
            ),
            exception=RDFError,
            function=RDF,
            traj=TrajectoryReader("traj.xyz"),
            reference_species=["Na"],
            target_species=["X"],
            delta_r=0.1,
        )

    def test__init__delta_r_zero(self, caplog):
        system = AtomicSystem(
            atoms=[Atom("h")],
            pos=np.array([[0, 0, 0]]),
            cell=Cell(10, 10, 10, 90, 90, 90)
        )
        traj = Trajectory([system])

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
            "The delta_r value of the RDF analysis has to be "
            "greater than zero - it actually is 0.0!"
            ),
            exception=RDFError,
            function=RDF,
            traj=traj,
            reference_species=["h"],
            target_species=["h"],
            delta_r=0.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
            "The delta_r value of the RDF analysis has to be "
            "greater than zero - it actually is 0.0!"
            ),
            exception=RDFError,
            function=RDF,
            traj=traj,
            reference_species=["h"],
            target_species=["h"],
            n_bins=5,
            delta_r=0.0,
        )

    def test_run_vacuum_trajectory(self, caplog):
        system1 = AtomicSystem(
            atoms=[Atom("h"), Atom("h")],
            pos=np.array([[0, 0, 0], [1, 0, 0]]),
            cell=Cell()
        )
        system2 = AtomicSystem(
            atoms=[Atom("h"), Atom("h")],
            pos=np.array([[0, 0, 0], [1, 0, 0]]),
            cell=Cell()
        )
        traj = Trajectory([system1, system2])

        rdf = RDF(traj, ["h"], ["h"], delta_r=1.0, r_max=5.0)

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
            "The provided trajectory is in vacuum, so the "
            "normalization of the RDF analysis requires a "
            "finite cell volume. Please provide a trajectory "
            "with box information."
            ),
            exception=RDFError,
            function=rdf.run,
        )

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test__init__uses_first_frame_topology_without_reader_topology(
        self, test_with_data_dir
    ):
        traj_reader = TrajectoryReader("traj.xyz")

        rdf = RDF(traj_reader, ["X"], ["X"], delta_r=0.1, n_bins=5)

        assert rdf.topology == rdf.first_frame.topology
        assert rdf.reference_indices.tolist() == [0, 1]
        assert rdf.target_indices.tolist() == [0, 1]

        _bin_middle_points, normalized_bins, *_ = rdf.run()

        assert np.isfinite(normalized_bins).all()

    def test_run_with_no_intra_molecular(self):
        rdf = RDF(
            _make_no_intra_trajectory(),
            ["H"],
            ["H"],
            delta_r=0.5,
            n_bins=5,
            no_intra_molecular=True
        )

        (
            _bin_middle_points,
            normalized_bins,
            integrated_bins,
            normalized_bins2,
            differential_bins
        ) = rdf.run()

        assert len(rdf.target_index_combinations) == 2
        assert np.array_equal(rdf.target_index_combinations[0], np.array([1]))
        assert np.array_equal(rdf.target_index_combinations[1], np.array([0]))
        assert np.allclose(integrated_bins, np.array([0.0, 0.0, 1.0, 1.0, 1.0]))
        assert np.isfinite(normalized_bins).all()
        assert np.isfinite(normalized_bins2).all()
        assert np.isfinite(differential_bins).all()

    def test_no_intra_molecular_excludes_own_residue(self):
        # two water molecules (atoms 0-5) followed by eight sodium
        # residues, so that every reference index is smaller than the
        # number of residues and no IndexError can hide the bug.
        rdf = RDF(
            _make_residue_trajectory(n_water=2, n_sodium=8),
            ["O"],
            ["H"],
            delta_r=0.5,
            n_bins=12,
            no_intra_molecular=True
        )

        assert rdf.reference_indices.tolist() == [0, 3]
        assert rdf.target_indices.tolist() == [1, 2, 4, 5]

        rdf.run()

        # each oxygen has to lose exactly the hydrogens of its own water
        assert len(rdf.target_index_combinations) == 2
        assert rdf.target_index_combinations[0].tolist() == [4, 5]
        assert rdf.target_index_combinations[1].tolist() == [1, 2]

        # hand-computed inter-molecular distances:
        # O(0) - H(4) = 5.3          -> bin 10
        # O(0) - H(5) = sqrt(19.49)  -> bin 8
        # O(3) - H(1) = 3.3          -> bin 6
        # O(3) - H(2) = sqrt(19.49)  -> bin 8
        # the intra-molecular O-H distances of 1.0 (bin 2) are excluded
        expected_bins = np.array([0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0])
        assert np.array_equal(rdf.bins, expected_bins)

    def test_no_intra_molecular_with_more_atoms_than_residues(self):
        # three water molecules, i.e. nine atoms but only three residues
        rdf = RDF(
            _make_residue_trajectory(n_water=3, n_sodium=0),
            ["O"],
            ["H"],
            delta_r=0.5,
            n_bins=12,
            no_intra_molecular=True
        )

        assert rdf.reference_indices.tolist() == [0, 3, 6]

        rdf.run()

        assert rdf.target_index_combinations[0].tolist() == [4, 5, 7, 8]
        assert rdf.target_index_combinations[1].tolist() == [1, 2, 7, 8]
        assert rdf.target_index_combinations[2].tolist() == [1, 2, 4, 5]

    def test_run_skips_self_pairs_for_overlapping_selections(self):
        system = AtomicSystem(
            atoms=[Atom("H"), Atom("H")],
            pos=np.array([[0, 0, 0], [1, 0, 0]]),
            cell=Cell(10, 10, 10, 90, 90, 90)
        )

        rdf = RDF(Trajectory([system]), ["H"], ["H"], delta_r=0.5, n_bins=4)

        (
            _bin_middle_points,
            normalized_bins,
            integrated_bins,
            normalized_bins2,
            differential_bins
        ) = rdf.run()

        np.testing.assert_allclose(rdf.bins, np.array([0.0, 0.0, 2.0, 0.0]))
        np.testing.assert_allclose(
            integrated_bins,
            np.array([0.0, 0.0, 1.0, 1.0])
        )
        assert np.isfinite(normalized_bins).all()
        assert np.isfinite(normalized_bins2).all()
        assert np.isfinite(differential_bins).all()

    def test_matches_ase_partial_rdf_reference(self):
        from ase import Atoms
        from ase.geometry.rdf import get_rdf as ase_get_rdf

        (
            symbols,
            box_length,
            positions_by_frame,
            trajectory,
        ) = _make_partial_rdf_reference_trajectory()
        delta_r = 0.5
        r_max = 5.0
        n_bins = int(r_max / delta_r)

        bin_centers, normalized_bins, *_ = RDF(
            trajectory,
            "H",
            "O",
            delta_r=delta_r,
            r_max=r_max,
        ).run()

        expected_centers = np.arange(delta_r / 2, r_max, delta_r)
        np.testing.assert_allclose(bin_centers, expected_centers)

        edges = np.linspace(0.0, r_max, n_bins + 1)
        h_indices = [i for i, symbol in enumerate(symbols) if symbol == "H"]
        o_indices = [i for i, symbol in enumerate(symbols) if symbol == "O"]
        counts = np.zeros(n_bins, dtype=float)
        ase_partial_bins = []

        for positions in positions_by_frame:
            atoms = Atoms(
                symbols=symbols,
                positions=positions,
                cell=[box_length, box_length, box_length],
                pbc=True,
            )
            distances = atoms.get_all_distances(mic=True)[
                np.ix_(h_indices, o_indices)
            ].ravel()
            distance_bins = np.floor_divide(distances, delta_r).astype(int)
            distance_bins = distance_bins[distances < r_max]
            counts += np.bincount(distance_bins, minlength=n_bins)[:n_bins]

            ase_rdf, _ase_distances = ase_get_rdf(
                atoms,
                rmax=r_max,
                nbins=n_bins,
                elements=[1, 8],
                no_dists=False,
            )
            ase_partial_bins.append(ase_rdf)

        shell_volumes = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
        target_density = len(o_indices) / box_length**3
        expected_bins = counts / (
            shell_volumes * target_density * len(h_indices) * len(positions_by_frame)
        )

        np.testing.assert_allclose(normalized_bins, expected_bins)
        ase_reference_bins = np.mean(ase_partial_bins, axis=0)
        # ASE changed partial RDF normalization across releases.
        if not np.allclose(ase_reference_bins, expected_bins):
            ase_reference_bins *= len(symbols) / len(o_indices)

        np.testing.assert_allclose(
            ase_reference_bins,
            expected_bins,
        )
        np.testing.assert_allclose(normalized_bins, ase_reference_bins)

    def test_run_twice_raises(self, caplog):
        """
        An RDF object is single-use: the first run returns the
        correct g(r), a second run raises a clear RDFError before
        the accumulators are touched again (previously the second
        call silently accumulated into self.bins and returned a
        doubled g(r)).
        """
        rdf = RDF(
            _make_no_intra_trajectory(), ["H"], ["H"], delta_r=0.5, n_bins=5
        )
        _, normalized_bins, *_ = rdf.run()

        reference = RDF(
            _make_no_intra_trajectory(), ["H"], ["H"], delta_r=0.5, n_bins=5
        )
        _, reference_bins, *_ = reference.run()

        assert np.allclose(normalized_bins, reference_bins)

        bins_after_first_run = rdf.bins.copy()

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="RDF",
            logging_level="ERROR",
            message_to_test=(
                "This RDF analysis object has already been run; "
                "construct a new one to run the analysis again."
            ),
            exception=RDFError,
            function=rdf.run,
        )

        # the guard preempts the accumulation, so the bin counts
        # are not doubled by the failed second call
        assert np.array_equal(rdf.bins, bins_after_first_run)
