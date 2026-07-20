"""
A module to test the ADF class (setup, angle calculation and output).
"""

import numpy as np
import pytest

from PQAnalysis.analysis.adf.exceptions import ADFError
from PQAnalysis.analysis import ADF
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



def _triangle_trajectory(pos):
    """One frame with a center O and two ligand H atoms."""
    return Trajectory([
        AtomicSystem(
            atoms=[Atom("O"), Atom("H"), Atom("H")],
            pos=np.array(pos, dtype=np.float64),
            cell=Cell(20.0, 20.0, 20.0),
        )
    ])


def _run_counts(traj, *args, **kwargs):
    adf = ADF(traj, *args, **kwargs)
    _mids, _norm, counts, _sinc = adf.run()
    return adf, counts



class TestADFSetup:

    def test_default_n_angle_bins(self):
        traj = _triangle_trajectory([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        adf = ADF(traj, "O", "H")

        assert adf.n_angle_bins == 180
        assert np.isclose(adf.delta_angle, 1.0)
        assert np.isclose(adf.max_angle, 180.0)
        assert len(adf.angle_bin_middle_points) == 180
        assert np.isclose(adf.angle_bin_middle_points[0], 0.5)
        assert np.isclose(adf.angle_bin_middle_points[-1], 179.5)

    def test_delta_angle_only(self):
        traj = _triangle_trajectory([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        adf = ADF(traj, "O", "H", delta_angle=2.0)

        assert adf.n_angle_bins == 90
        assert np.isclose(adf.delta_angle, 2.0)

    def test_n_angle_bins_only(self):
        traj = _triangle_trajectory([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        adf = ADF(traj, "O", "H", n_angle_bins=360)

        assert adf.n_angle_bins == 360
        assert np.isclose(adf.delta_angle, 0.5)

    def test_both_bins_raises(self, caplog):
        traj = _triangle_trajectory([[0, 0, 0], [1, 0, 0], [0, 1, 0]])

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ADF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
            "It is not possible to specify both n_angle_bins and "
            "delta_angle in the same ADF analysis as this would "
            "lead to ambiguous results."
            ),
            exception=ADFError,
            function=ADF,
            traj=traj,
            reference_species="O",
            target_species="H",
            n_angle_bins=180,
            delta_angle=1.0,
        )

    def test_delta_angle_too_large_raises(self, caplog):
        traj = _triangle_trajectory([[0, 0, 0], [1, 0, 0], [0, 1, 0]])

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ADF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
            "The provided delta_angle is larger than the angle range "
            "of 180 degrees, resulting in no angle bins."
            ),
            exception=ADFError,
            function=ADF,
            traj=traj,
            reference_species="O",
            target_species="H",
            delta_angle=200.0,
        )

    def test_empty_trajectory_raises(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ADF.__qualname__,
            logging_level="ERROR",
            message_to_test="Trajectory cannot be of length 0.",
            exception=ADFError,
            function=ADF,
            traj=Trajectory(),
            reference_species="O",
            target_species="H",
        )

    def test_target_selection_2_defaults_to_target(self):
        traj = _triangle_trajectory([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        adf = ADF(traj, "O", "H")

        assert adf.target_species_2 == "H"
        assert np.array_equal(adf.target_indices, adf.target_indices_2)

    def test_n_frames_and_n_atoms(self):
        traj = _triangle_trajectory([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        adf = ADF(traj, "O", "H")

        assert adf.n_frames == 1
        assert adf.n_atoms == 3

    def test_setup_bin_middle_points(self):
        mids = ADF._setup_bin_middle_points(4, 45.0)
        assert np.allclose(mids, np.array([22.5, 67.5, 112.5, 157.5]))



class TestADFAngles:

    def test_right_angle(self):
        # O at origin, H on +x and +y -> j-i-k angle is 90 degrees
        traj = _triangle_trajectory([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        adf, counts = _run_counts(traj, "O", "H", n_angle_bins=180)

        assert counts.sum() == 2  # ordered pairs (j,k) and (k,j)
        assert counts[90] == 2
        assert np.count_nonzero(counts) == 1

    def test_tetrahedral_angle(self):
        a = 1.0 / np.sqrt(3.0)
        traj = _triangle_trajectory(
            [[0, 0, 0], [a, a, a], [a, -a, -a]]
        )
        adf, counts = _run_counts(traj, "O", "H", n_angle_bins=180)

        # 109.47 degrees -> bin 109 (covering [109, 110))
        assert counts[109] == 2
        assert np.count_nonzero(counts) == 1

    def test_straight_angle_is_dropped(self):
        # exact 180 degrees lands in bin n_bins and is discarded
        traj = _triangle_trajectory([[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
        adf, counts = _run_counts(traj, "O", "H", n_angle_bins=180)

        assert counts.sum() == 0

    def test_ordered_pairs_counted_once_with_distinct_ligand_sets(self):
        # C center, distinct ligand-1 (H) and ligand-2 (O) sets -> the
        # angle is counted once per (j, k) pair, not doubled
        traj = Trajectory([
            AtomicSystem(
                atoms=[Atom("C"), Atom("H"), Atom("O")],
                pos=np.array([[0, 0, 0], [1.0, 0, 0], [0, 1.0, 0]]),
                cell=Cell(20.0, 20.0, 20.0),
            )
        ])
        adf, counts = _run_counts(
            traj, "C", "H", target_species_2="O", n_angle_bins=180
        )

        assert counts.sum() == 1
        assert counts[90] == 1

    def test_radial_gate_excludes_far_ligands(self):
        # one near H (r = 1) and one far H (r = 5) around the center;
        # gating i-j to [0, 2) keeps only the near ligand as j, so no
        # (j, k) pair with j != k survives within the near shell
        traj = Trajectory([
            AtomicSystem(
                atoms=[Atom("O"), Atom("H"), Atom("H"), Atom("H")],
                pos=np.array(
                    [[0, 0, 0], [1.0, 0, 0], [0, 1.0, 0], [5.0, 0, 0]]
                ),
                cell=Cell(20.0, 20.0, 20.0),
            )
        ])

        _adf_all, counts_all = _run_counts(traj, "O", "H", n_angle_bins=180)
        _adf_gate, counts_gate = _run_counts(
            traj, "O", "H", n_angle_bins=180, r_max_1=2.0, r_max_2=2.0
        )

        # without a gate the far ligand contributes; with the gate the
        # total number of counted triplets strictly decreases
        assert counts_gate.sum() < counts_all.sum()
        assert counts_gate.sum() > 0

    def test_normalized_integral_is_one(self):
        rng = np.random.default_rng(7)
        pos = rng.uniform(0.0, 12.0, size=(20, 3))
        atoms = [Atom("O") if i % 2 == 0 else Atom("H") for i in range(20)]
        traj = Trajectory([
            AtomicSystem(atoms=atoms, pos=pos, cell=Cell(12.0, 12.0, 12.0))
        ])

        adf = ADF(traj, "O", "H", n_angle_bins=180)
        _mids, norm, counts, sinc = adf.run()

        assert np.isclose((norm * adf.delta_angle).sum(), 1.0)
        assert np.isclose((sinc * adf.delta_angle).sum(), 1.0)
        assert np.isclose(counts.sum(), adf.angle_bins.sum())

    def test_sin_correction_flattens_isotropic_distribution(self):
        # isotropic ligand directions -> the raw ADF is biased towards
        # 90 degrees (the sin(theta) solid-angle factor); the
        # sin-corrected density must be markedly flatter
        rng = np.random.default_rng(2026)
        n_ligands = 150

        directions = rng.normal(size=(n_ligands, 3))
        directions /= np.linalg.norm(directions, axis=1, keepdims=True)
        radii = rng.uniform(1.0, 3.0, size=(n_ligands, 1))

        center = np.array([20.0, 20.0, 20.0])
        ligand_positions = center + directions * radii

        pos = np.vstack([center, ligand_positions])
        atoms = [Atom("O")] + [Atom("H")] * n_ligands

        traj = Trajectory([
            AtomicSystem(atoms=atoms, pos=pos, cell=Cell(40.0, 40.0, 40.0))
        ])

        adf = ADF(traj, "O", "H", n_angle_bins=18)
        _mids, norm, _counts, sinc = adf.run()

        # interior bins (avoid the sparsely populated extreme angles)
        interior = slice(2, 16)

        def coefficient_of_variation(values):
            return np.std(values) / np.mean(values)

        assert (
            coefficient_of_variation(sinc[interior])
            < coefficient_of_variation(norm[interior])
        )



class TestADFNormalizationHelpers:

    def test_normalize_empty_counts(self):
        counts = np.zeros(5)
        assert np.array_equal(ADF._normalize(counts, 1.0), np.zeros(5))

    def test_normalize_integral(self):
        counts = np.array([1.0, 3.0, 6.0])
        delta = 2.0
        normalized = ADF._normalize(counts, delta)

        assert np.isclose((normalized * delta).sum(), 1.0)

    def test_sin_correct_empty_counts(self):
        counts = np.zeros(5)
        mids = ADF._setup_bin_middle_points(5, 36.0)
        assert np.array_equal(
            ADF._sin_correct(counts, mids, 36.0), np.zeros(5)
        )

    def test_sin_correct_integral(self):
        counts = np.array([2.0, 5.0, 5.0, 2.0])
        mids = ADF._setup_bin_middle_points(4, 45.0)
        delta = 45.0
        corrected = ADF._sin_correct(counts, mids, delta)

        assert np.isclose((corrected * delta).sum(), 1.0)



class TestADFDispatch:

    @staticmethod
    def _write_traj(path, n_frames=4, n_atoms=15, box=12.0, seed=1):
        rng = np.random.default_rng(seed)
        with open(path, "w", encoding="utf-8") as file:
            for _ in range(n_frames):
                pos = rng.uniform(0.0, box, size=(n_atoms, 3))
                file.write(f"{n_atoms} {box} {box} {box}\n\n")
                for i, (x, y, z) in enumerate(pos):
                    name = "O" if i % 3 == 0 else "H"
                    file.write(f"{name} {x:.6f} {y:.6f} {z:.6f}\n")
        return str(path)

    def test_reader_uses_raw_fast_path(self, tmp_path):
        filename = self._write_traj(tmp_path / "traj.xyz")

        adf = ADF(TrajectoryReader(filename), "O", "H", n_angle_bins=180)
        assert adf._raw_reader is not None
        assert adf.frame_generator is None

    def test_in_memory_uses_original_path(self, tmp_path):
        filename = self._write_traj(tmp_path / "traj.xyz")

        adf = ADF(
            TrajectoryReader(filename).read(), "O", "H", n_angle_bins=180
        )
        assert adf._raw_reader is None
        assert adf.frame_generator is not None

    def test_fast_path_matches_original_path(self, tmp_path):
        filename = self._write_traj(tmp_path / "traj.xyz")

        adf_fast = ADF(TrajectoryReader(filename), "O", "H", n_angle_bins=180)
        result_fast = np.column_stack(adf_fast.run())

        adf_memory = ADF(
            TrajectoryReader(filename).read(), "O", "H", n_angle_bins=180
        )
        result_memory = np.column_stack(adf_memory.run())

        assert adf_fast._raw_reader is not None
        assert adf_memory._raw_reader is None
        assert np.array_equal(result_fast, result_memory)
        assert adf_fast.angle_bins.sum() > 0

    def test_fast_path_frame_with_missing_atoms(
        self, tmp_path, caplog, monkeypatch
    ):
        filename = self._write_traj(tmp_path / "traj.xyz")

        adf_fast = ADF(TrajectoryReader(filename), "O", "H", n_angle_bins=180)

        def _truncated_generator():
            yield (
                np.zeros((3, 3), dtype=np.float32),
                Cell(12.0, 12.0, 12.0),
            )

        monkeypatch.setattr(
            adf_fast._raw_reader,
            "raw_frame_generator",
            _truncated_generator,
        )

        max_index = int(
            max(
                adf_fast.reference_indices.max(),
                adf_fast.target_indices.max(),
                adf_fast.target_indices_2.max(),
            )
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ADF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
                "Frame 1 of the trajectory provides only 3 atoms, "
                f"but the selections reference the atom index {max_index}. "
                "Please provide a trajectory with a consistent "
                "number of atoms."
            ),
            exception=ADFError,
            function=adf_fast.run,
        )



class TestADFTypeChecking:

    def test_init_type_checking(self, caplog):
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
            function=ADF,
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
            function=ADF,
            traj=Trajectory(),
            reference_species=Trajectory(),
            target_species=["h"],
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "n_angle_bins",
            -1,
            PositiveInt | None
            ),
            exception=PQTypeError,
            function=ADF,
            traj=Trajectory(),
            reference_species=["h"],
            target_species=["h"],
            n_angle_bins=-1,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "delta_angle",
            -1,
            PositiveReal | None
            ),
            exception=PQTypeError,
            function=ADF,
            traj=Trajectory(),
            reference_species=["h"],
            target_species=["h"],
            delta_angle=-1,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "r_min_1",
            -1,
            PositiveReal | None
            ),
            exception=PQTypeError,
            function=ADF,
            traj=Trajectory(),
            reference_species=["h"],
            target_species=["h"],
            r_min_1=-1,
        )
