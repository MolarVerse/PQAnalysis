"""
Tests for the MSD class.
"""

import sys

import numpy as np
import pytest

from PQAnalysis import config
from PQAnalysis.analysis.msd import MSD
from PQAnalysis.analysis.msd.exceptions import MSDError
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.traj import Trajectory, TrajectoryFormat

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging, assert_logging_with_exception

#: The module defining the MSD class (the package attribute ``msd``
#: is shadowed by the api function of the same name).
msd_module = sys.modules[MSD.__module__]

# pylint: disable=protected-access



def _make_trajectory(positions, box=100.0, names=None):
    """
    Builds a trajectory from an (n_frames, n_atoms, 3) position array.
    """
    positions = np.asarray(positions, dtype=float)
    n_atoms = positions.shape[1]

    if names is None:
        names = ["O"] * n_atoms

    cell = Cell(box, box, box, 90, 90, 90)
    atoms = [Atom(name) for name in names]

    systems = [
        AtomicSystem(atoms=atoms, pos=frame_positions, cell=cell)
        for frame_positions in positions
    ]

    return Trajectory(systems)



def _reference_msd(wrapped, window, gap, box, n_start=0):
    """
    Brute-force emulation of the legacy Diffcalc algorithm.
    """
    n_frames, n_atoms, _ = wrapped.shape

    stop_frame = (n_frames - window) // gap * gap
    total_origins = stop_frame // gap

    # origins spawn at multiples of the gap not smaller than n_start
    first_origin = gap * -(-max(n_start, 1) // gap)

    msd = np.zeros((window + 1, 3))

    for origin in range(first_origin, stop_frame + 1, gap):
        image = np.zeros((n_atoms, 3))

        for frame in range(origin, origin + window + 1):
            if frame > origin:
                displacement = wrapped[frame - 1] - wrapped[frame - 2]
                image -= box * np.rint(displacement / box)

            total = wrapped[frame - 1] - wrapped[origin - 1] + image
            msd[frame - origin] += (total**2).sum(axis=0)

    return msd / (n_atoms * total_origins)



class TestMSD:

    """
    Tests for the MSD class.
    """

    def test_two_frame_msd(self):
        positions = [
            [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
            [[0.5, 0.0, 0.0], [1.0, 1.0, 2.0]],
        ]
        traj = _make_trajectory(positions)

        msd = MSD(traj, "O", window=1, gap=1)
        lags, msd_x, msd_y, msd_z, msd_tot = msd.run()

        assert np.array_equal(lags, [0, 1])
        assert np.allclose(msd_x, [0.0, 0.125])
        assert np.allclose(msd_y, [0.0, 0.0])
        assert np.allclose(msd_z, [0.0, 0.5])
        assert np.allclose(msd_tot, [0.0, 0.625])

    def test_three_frame_msd(self):
        positions = [
            [[0.0, 0.0, 0.0]],
            [[1.0, 0.5, 0.0]],
            [[3.0, 1.5, 0.0]],
        ]
        traj = _make_trajectory(positions)

        msd = MSD(traj, "O", window=2, gap=1)
        _, msd_x, msd_y, msd_z, _ = msd.run()

        # a single origin at the first frame contributes lags 0, 1 and 2
        assert np.allclose(msd_x, [0.0, 1.0, 9.0])
        assert np.allclose(msd_y, [0.0, 0.25, 2.25])
        assert np.allclose(msd_z, [0.0, 0.0, 0.0])

    def test_gap_one_full_window_matches_legacy(self, caplog):
        # legacy Diffcalc boundary case: gap == 1 with a trajectory
        # of exactly window frames spawns a single time origin at
        # the first frame instead of raising the
        # too-short-trajectory error; a warning is emitted instead
        positions = [
            [[0.0, 0.0, 0.0]],
            [[1.0, 0.5, 0.0]],
            [[3.0, 1.5, 0.0]],
        ]
        traj = _make_trajectory(positions)

        msd = assert_logging(
            caplog=caplog,
            logging_name="MSD",
            logging_level="WARNING",
            message_to_test=(
                "The trajectory contains exactly window = 3 frames "
                "with a gap of 1. Following the legacy Diffcalc "
                "convention a single time origin is spawned at the "
                "first frame. The final lag bin (lag = 3) can never "
                "be sampled and is written as exactly 0.0; a "
                "diffusion fit including this bin would be biased."
            ),
            function=MSD,
            traj=traj,
            target_species="O",
            window=3,
            gap=1,
        )

        assert msd.stop_frame == 1
        assert msd.total_origins == 1

        lags, msd_x, msd_y, msd_z, _ = msd.run()

        # pinned against the recompiled legacy Diffcalc binary
        # (window = 3; gap = 1; on the same three-frame trajectory):
        # the single origin fills lags 0 to 2, the final lag bin
        # is written as exactly 0.0
        assert np.array_equal(lags, [0, 1, 2, 3])
        assert np.allclose(msd_x, [0.0, 1.0, 9.0, 0.0])
        assert np.allclose(msd_y, [0.0, 0.25, 2.25, 0.0])
        assert np.allclose(msd_z, [0.0, 0.0, 0.0, 0.0])
        assert msd.msd_tot[3] == 0.0

    def test_unwrap_across_boundary(self):
        positions = [
            [[9.8, 0.1, 5.0]],
            [[0.1, 9.8, 5.0]],
        ]
        traj = _make_trajectory(positions, box=10.0)

        msd = MSD(traj, "O", window=1, gap=1)
        _, msd_x, msd_y, msd_z, _ = msd.run()

        # x crosses the upper boundary (+0.3), y the lower one (-0.3)
        assert np.allclose(msd_x, [0.0, 0.09])
        assert np.allclose(msd_y, [0.0, 0.09])
        assert np.allclose(msd_z, [0.0, 0.0])

    def test_matches_brute_force_reference(self):
        rng = np.random.default_rng(4711)
        n_atoms, n_frames, box = 6, 57, 8.0

        start = rng.uniform(0.0, box, size=(1, n_atoms, 3))
        steps = rng.normal(0.0, 0.4, size=(n_frames - 1, n_atoms, 3))
        unwrapped = np.concatenate(
            [start, start + np.cumsum(steps, axis=0)]
        )
        wrapped = unwrapped % box

        window, gap = 20, 5

        traj = _make_trajectory(wrapped, box=box)
        msd = MSD(traj, "O", window=window, gap=gap)
        _, msd_x, msd_y, msd_z, _ = msd.run()

        reference = _reference_msd(wrapped, window, gap, box)

        assert np.allclose(msd_x, reference[:, 0], rtol=1e-10, atol=1e-12)
        assert np.allclose(msd_y, reference[:, 1], rtol=1e-10, atol=1e-12)
        assert np.allclose(msd_z, reference[:, 2], rtol=1e-10, atol=1e-12)

    def test_selection_subset(self):
        rng = np.random.default_rng(1234)
        wrapped = rng.uniform(0.0, 5.0, size=(8, 4, 3))
        names = ["O", "H", "O", "H"]

        traj = _make_trajectory(wrapped, box=5.0, names=names)
        msd = MSD(traj, "O", window=2, gap=1)

        assert np.array_equal(msd.target_indices, [0, 2])

        _, msd_x, _, _, _ = msd.run()

        reference = _reference_msd(wrapped[:, [0, 2]], 2, 1, 5.0)

        assert np.allclose(msd_x, reference[:, 0], rtol=1e-10, atol=1e-12)

    def test_n_start_skips_leading_frames(self):
        rng = np.random.default_rng(99)
        wrapped = rng.uniform(0.0, 5.0, size=(30, 2, 3)) % 5.0

        traj = _make_trajectory(wrapped, box=5.0)
        msd = MSD(traj, "O", window=10, gap=5, n_start=6)

        _, msd_x, msd_y, msd_z, _ = msd.run()

        # origins may only spawn at multiples of the gap >= n_start
        # (here 10, 15 and 20), the legacy normalization by
        # total_origins = stop_frame // gap = 4 is kept nevertheless
        reference = _reference_msd(wrapped, 10, 5, 5.0, n_start=6)

        assert np.allclose(msd_x, reference[:, 0], rtol=1e-10, atol=1e-12)
        assert np.allclose(msd_y, reference[:, 1], rtol=1e-10, atol=1e-12)
        assert np.allclose(msd_z, reference[:, 2], rtol=1e-10, atol=1e-12)

    def test_run_with_time_step_sets_fit_results(self):
        rng = np.random.default_rng(7)
        wrapped = np.cumsum(
            rng.normal(0.0, 0.2, size=(40, 3, 3)), axis=0
        ) % 6.0

        traj = _make_trajectory(wrapped, box=6.0)
        msd = MSD(traj, "O", window=10, gap=5, time_step=0.5)

        assert msd.fit_window == 2

        msd.run()

        assert msd.fit_results is not None
        assert set(msd.fit_results.keys()) == {"x", "y", "z", "total"}

        for fit in msd.fit_results.values():
            assert np.isfinite(fit.diffusion_coefficient)
            assert np.isfinite(fit.slope)

    def test_fit_diffusion_linear(self):
        lags = np.arange(101)
        time_step = 0.1
        times = lags * time_step

        # per axis: MSD_axis = 2 * D_axis * t with slopes 0.6, 0.8, 1.0
        msd_x = 0.6 * times + 1.0
        msd_y = 0.8 * times + 2.0
        msd_z = 1.0 * times + 3.0
        msd_tot = msd_x + msd_y + msd_z

        results = MSD._fit_diffusion(
            lags,
            msd_x,
            msd_y,
            msd_z,
            msd_tot,
            time_step,
            fit_window=20
        )

        assert results["x"].slope == pytest.approx(0.6)
        assert results["y"].slope == pytest.approx(0.8)
        assert results["z"].slope == pytest.approx(1.0)
        assert results["total"].slope == pytest.approx(2.4)

        # D = slope / (2 * dim) * 1e-8 m^2/s
        assert results["x"].diffusion_coefficient == pytest.approx(0.3e-8)
        assert results["z"].diffusion_coefficient == pytest.approx(0.5e-8)
        assert results["total"].diffusion_coefficient == pytest.approx(0.4e-8)

        for fit in results.values():
            assert fit.r_squared == pytest.approx(1.0)
            assert fit.slope_stderr == pytest.approx(0.0, abs=1e-8)

    def test_fit_diffusion_uses_trailing_window(self):
        # piecewise data: the leading part (indices 0 to 80) has a
        # steep slope, only the exact trailing fit_window = 20 points
        # (indices 81 to 100) lie on the flat tail lines. Fitting the
        # leading window, ignoring fit_window (fitting everything) or
        # an off-by-one slice including index 80 (which is far off
        # the tail lines) all yield different slopes and r^2 < 1.
        lags = np.arange(101)
        time_step = 0.1
        times = lags * time_step
        fit_window = 20

        tail = lags > 80

        msd_x = np.where(tail, 0.6 * times + 100.0, 2.0 * times)
        msd_y = np.where(tail, 0.8 * times + 200.0, 3.0 * times)
        msd_z = np.where(tail, 1.0 * times + 300.0, 4.0 * times)
        msd_tot = msd_x + msd_y + msd_z

        results = MSD._fit_diffusion(
            lags,
            msd_x,
            msd_y,
            msd_z,
            msd_tot,
            time_step,
            fit_window=fit_window
        )

        assert results["x"].slope == pytest.approx(0.6)
        assert results["y"].slope == pytest.approx(0.8)
        assert results["z"].slope == pytest.approx(1.0)
        assert results["total"].slope == pytest.approx(2.4)

        assert results["x"].intercept == pytest.approx(100.0)
        assert results["y"].intercept == pytest.approx(200.0)
        assert results["z"].intercept == pytest.approx(300.0)
        assert results["total"].intercept == pytest.approx(600.0)

        assert results["x"].diffusion_coefficient == pytest.approx(0.3e-8)
        assert results["total"].diffusion_coefficient == pytest.approx(0.4e-8)

        # r^2 == 1 and stderr == 0 hold only for the exact trailing
        # window: including even one leading point breaks both (an
        # off-by-one slice gives stderr of order 1 and r^2 well
        # below 1; the tolerance only allows for float rounding of
        # the large intercepts)
        for fit in results.values():
            assert fit.r_squared == pytest.approx(1.0, abs=1e-9)
            assert fit.slope_stderr == pytest.approx(0.0, abs=1e-6)

    def test_fit_diffusion_noisy(self):
        rng = np.random.default_rng(42)

        lags = np.arange(101)
        time_step = 0.1
        times = lags * time_step

        noise = rng.normal(0.0, 0.05, size=(4, lags.size))

        msd_x = 1.2 * times + noise[0]
        msd_y = 1.2 * times + noise[1]
        msd_z = 1.2 * times + noise[2]
        msd_tot = 3.6 * times + noise[3]

        results = MSD._fit_diffusion(
            lags,
            msd_x,
            msd_y,
            msd_z,
            msd_tot,
            time_step,
            fit_window=30
        )

        for fit in results.values():
            assert 0.0 < fit.r_squared < 1.0
            assert fit.slope_stderr > 0.0
            assert fit.diffusion_coefficient_stderr > 0.0

        assert results["x"].slope == pytest.approx(1.2, rel=0.1)
        assert results["total"].slope == pytest.approx(3.6, rel=0.1)

    def test_window_not_multiple_of_gap(self, caplog):
        traj = _make_trajectory(np.zeros((30, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSD",
            logging_level="ERROR",
            message_to_test=(
                "The window size 15 has to be an integer multiple of the gap 4."
            ),
            exception=MSDError,
            function=MSD,
            traj=traj,
            target_species="O",
            window=15,
            gap=4,
        )

    def test_trajectory_too_short(self, caplog):
        traj = _make_trajectory(np.zeros((12, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSD",
            logging_level="ERROR",
            message_to_test=(
                "The trajectory with 12 frames is too short to establish "
                "a window of 10 frames with a gap of 5 frames. At least "
                "window + gap = 15 frames are required (or exactly "
                "window frames for gap == 1)."
            ),
            exception=MSDError,
            function=MSD,
            traj=traj,
            target_species="O",
            window=10,
            gap=5,
        )

    def test_trajectory_too_short_gap_one(self, caplog):
        # for gap == 1 only a trajectory of exactly window frames is
        # accepted as the legacy single-origin boundary case - any
        # shorter trajectory has to raise
        traj = _make_trajectory(np.zeros((8, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSD",
            logging_level="ERROR",
            message_to_test=(
                "The trajectory with 8 frames is too short to establish "
                "a window of 10 frames with a gap of 1 frames. At least "
                "window + gap = 11 frames are required (or exactly "
                "window frames for gap == 1)."
            ),
            exception=MSDError,
            function=MSD,
            traj=traj,
            target_species="O",
            window=10,
            gap=1,
        )

    def test_empty_trajectory(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSD",
            logging_level="ERROR",
            message_to_test="Trajectory cannot be of length 0.",
            exception=MSDError,
            function=MSD,
            traj=Trajectory(),
            target_species="O",
            window=1,
            gap=1,
        )

    def test_negative_n_start(self, caplog):
        traj = _make_trajectory(np.zeros((30, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSD",
            logging_level="ERROR",
            message_to_test="n_start must be a non-negative integer.",
            exception=MSDError,
            function=MSD,
            traj=traj,
            target_species="O",
            window=10,
            gap=5,
            n_start=-1,
        )

    def test_n_start_too_large(self, caplog):
        traj = _make_trajectory(np.zeros((30, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSD",
            logging_level="ERROR",
            message_to_test=(
                "The starting frame 25 is too large: time origins only "
                "spawn at multiples of the gap 5 up to stop_frame = "
                "(n_frames - window) // gap * gap = 20 (n_frames = 30, "
                "window = 10), so no time origin could spawn."
            ),
            exception=MSDError,
            function=MSD,
            traj=traj,
            target_species="O",
            window=10,
            gap=5,
            n_start=25,
        )

    def test_n_start_between_stop_frame_and_window_limit(self, caplog):
        # n_start = 5 satisfies n_start <= n_frames - window = 6 but
        # exceeds stop_frame = (14 - 8) // 4 * 4 = 4, so no origin
        # could ever spawn - this used to run through silently and
        # write an all-zero MSD
        traj = _make_trajectory(np.zeros((14, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSD",
            logging_level="ERROR",
            message_to_test=(
                "The starting frame 5 is too large: time origins only "
                "spawn at multiples of the gap 4 up to stop_frame = "
                "(n_frames - window) // gap * gap = 4 (n_frames = 14, "
                "window = 8), so no time origin could spawn."
            ),
            exception=MSDError,
            function=MSD,
            traj=traj,
            target_species="O",
            window=8,
            gap=4,
            n_start=5,
        )

    def test_empty_selection(self, caplog):
        traj = _make_trajectory(np.zeros((30, 1, 3)), names=["H"])

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSD",
            logging_level="ERROR",
            message_to_test="The target selection does not select any atoms.",
            exception=MSDError,
            function=MSD,
            traj=traj,
            target_species="O",
            window=10,
            gap=5,
        )

    def test_zero_time_step(self, caplog):
        # time_step = 0.0 passes the PositiveReal type check but
        # would make every fit abscissa identical - it has to be
        # rejected before any frame is streamed
        traj = _make_trajectory(np.zeros((30, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSD",
            logging_level="ERROR",
            message_to_test=(
                "The time_step must be a positive real number."
            ),
            exception=MSDError,
            function=MSD,
            traj=traj,
            target_species="O",
            window=10,
            gap=5,
            time_step=0.0,
        )

    def test_fit_window_too_small(self, caplog):
        # a linear regression needs at least two points, otherwise
        # scipy returns an all-NaN fit
        traj = _make_trajectory(np.zeros((30, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSD",
            logging_level="ERROR",
            message_to_test=(
                "The fit_window 1 must be at least 2 to perform a "
                "linear diffusion fit."
            ),
            exception=MSDError,
            function=MSD,
            traj=traj,
            target_species="O",
            window=10,
            gap=5,
            fit_window=1,
        )

    def test_fit_window_too_large(self, caplog):
        traj = _make_trajectory(np.zeros((30, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSD",
            logging_level="ERROR",
            message_to_test=(
                "The fit_window 12 cannot be larger than window + 1 = 11."
            ),
            exception=MSDError,
            function=MSD,
            traj=traj,
            target_species="O",
            window=10,
            gap=5,
            fit_window=12,
        )

    def test_frame_atom_count_mismatch_multi_file(self, tmp_path, caplog):
        # multiple trajectory files with different atom counts used
        # to be silently minimum-imaged into garbage MSD values
        file_a = tmp_path / "a.xyz"
        file_b = tmp_path / "b.xyz"

        file_a.write_text(
            (
                "2 100.0 100.0 100.0\n\n"
                "O 0.0 0.0 0.0\nH 1.0 0.0 0.0\n"
                "2 100.0 100.0 100.0\n\n"
                "O 0.1 0.0 0.0\nH 1.0 0.0 0.0\n"
            ),
            encoding="utf-8",
        )
        file_b.write_text(
            (
                "1 100.0 100.0 100.0\n\n"
                "O 0.2 0.0 0.0\n"
                "1 100.0 100.0 100.0\n\n"
                "O 0.3 0.0 0.0\n"
            ),
            encoding="utf-8",
        )

        reader = TrajectoryReader([str(file_a), str(file_b)])
        msd = MSD(reader, "O", window=2, gap=1)

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSD",
            logging_level="ERROR",
            message_to_test=(
                "Frame 3 of the trajectory does not provide positions "
                "for all 2 atoms of the topology. Please provide a "
                "position trajectory (e.g. .xyz files) with a "
                "consistent number of atoms."
            ),
            exception=MSDError,
            function=msd.run,
        )

    def test_extxyz_reader_uses_non_raw_old_path(self, tmp_path):
        # a TrajectoryReader is only routed to the raw fast path for
        # plain XYZ trajectories; every other frame-carrying format
        # (here EXTXYZ) is streamed through the per-frame
        # AtomicSystem frame_generator (the old non-fast path). The
        # result has to be identical to running the in-memory
        # trajectory of the same wrapped coordinates and cell.
        rng = np.random.default_rng(2024)
        box = (10.0, 12.0, 14.0)
        n_frames, n_atoms = 8, 3

        wrapped = rng.uniform(0.0, 5.0, size=(n_frames, n_atoms, 3))

        file = tmp_path / "traj.extxyz"
        lines = []
        for frame in wrapped:
            lines.append(str(n_atoms))
            lines.append(
                f'Lattice="{box[0]} 0 0 0 {box[1]} 0 0 0 {box[2]}" '
                "Properties=species:S:1:pos:R:3"
            )
            for atom in frame:
                lines.append(f"O {atom[0]} {atom[1]} {atom[2]}")
        file.write_text("\n".join(lines) + "\n", encoding="utf-8")

        reader = TrajectoryReader(str(file))

        assert reader.traj_format == TrajectoryFormat.EXTXYZ

        msd = MSD(reader, "O", window=2, gap=1)

        # the EXTXYZ reader is not eligible for the raw fast path:
        # it lands on the elif TrajectoryReader (old) branch that sets
        # up the per-frame frame_generator instead of a _raw_reader
        assert msd._raw_reader is None
        assert msd.frame_generator is not None
        assert msd.n_frames == n_frames

        _, msd_x, msd_y, msd_z, _ = msd.run()

        cell = Cell(*box, 90, 90, 90)
        atoms = [Atom("O")] * n_atoms
        systems = [
            AtomicSystem(atoms=atoms, pos=frame, cell=cell)
            for frame in wrapped
        ]
        reference = MSD(Trajectory(systems), "O", window=2, gap=1)
        _, ref_x, ref_y, ref_z, _ = reference.run()

        assert np.allclose(msd_x, ref_x, rtol=1e-5, atol=1e-6)
        assert np.allclose(msd_y, ref_y, rtol=1e-5, atol=1e-6)
        assert np.allclose(msd_z, ref_z, rtol=1e-5, atol=1e-6)

    def test_velocity_trajectory_rejected(self, caplog):
        # frames of a velocity trajectory carry no positions and
        # used to crash with a raw IndexError
        cell = Cell(100.0, 100.0, 100.0, 90, 90, 90)
        atoms = [Atom("O"), Atom("O")]

        systems = [
            AtomicSystem(atoms=atoms, vel=np.ones((2, 3)), cell=cell)
            for _ in range(3)
        ]
        traj = Trajectory(systems)

        msd = MSD(traj, "O", window=1, gap=1)

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSD",
            logging_level="ERROR",
            message_to_test=(
                "Frame 1 of the trajectory does not provide positions "
                "for all 2 atoms of the topology. Please provide a "
                "position trajectory (e.g. .xyz files) with a "
                "consistent number of atoms."
            ),
            exception=MSDError,
            function=msd.run,
        )

    def test_progress_bar_binds_config_at_call_time(self, monkeypatch):
        # config.with_progress_bar is set by the CLI after the module
        # import, so it must be read at call time, not bound by value
        # at import time
        captured = {}

        def fake_tqdm(iterable, **kwargs):
            captured.update(kwargs)
            return iterable

        monkeypatch.setattr(msd_module, "tqdm", fake_tqdm)

        monkeypatch.setattr(config, "with_progress_bar", False)
        MSD(_make_trajectory(np.zeros((4, 1, 3))), "O", window=2, gap=1).run()
        assert captured["disable"] is True

        captured.clear()

        monkeypatch.setattr(config, "with_progress_bar", True)
        MSD(_make_trajectory(np.zeros((4, 1, 3))), "O", window=2, gap=1).run()
        assert captured["disable"] is False

    def test_defaults(self):
        traj = _make_trajectory(np.zeros((1200, 1, 3)))

        msd = MSD(traj, "O")

        assert msd.window == 1000
        assert msd.gap == 10
        assert msd.n_start == 0
        assert msd.time_step is None
        assert msd.fit_window == 200
        assert msd.n_origins_max == 100
        assert msd.stop_frame == 200
        assert msd.total_origins == 20

    def test_unwrap_shift_vacuum(self):
        displacement = np.array([[1.0, 2.0, 3.0]])

        shift = MSD._unwrap_shift(displacement, Cell())

        assert np.array_equal(shift, np.zeros((1, 3)))

    def test_unwrap_shift_orthorhombic(self):
        cell = Cell(10.0, 12.0, 14.0, 90, 90, 90)
        displacement = np.array([[9.7, -11.9, 3.0]])

        shift = MSD._unwrap_shift(displacement, cell)

        assert np.allclose(shift, [[-10.0, 12.0, 0.0]])
