import numpy as np
import pytest

from PQAnalysis.analysis.rdf.exceptions import RDFError
from PQAnalysis.analysis import RDF
from PQAnalysis.traj import Trajectory
from PQAnalysis.core import Atom, Cell
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.topology import SelectionCompatible
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

        system1 = AtomicSystem(cell=Cell(10, 10, 10, 90, 90, 90))
        system2 = AtomicSystem(cell=Cell(16, 13, 12, 90, 90, 90))

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
