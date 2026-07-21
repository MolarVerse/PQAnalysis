"""
Tests for the GreenKubo class and its module-level helpers.
"""

import sys

import numpy as np
import pytest

from PQAnalysis import config
from PQAnalysis.analysis.green_kubo import (
    GreenKubo,
    ANGSTROM2_PER_S2_PS_TO_M2_PER_S,
    M2_PER_S_TO_CM2_PER_S,
    cumulative_trapezoid,
    velocity_acf_direct,
    velocity_acf_fft,
)
from PQAnalysis.analysis.green_kubo.exceptions import GreenKuboError
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.traj import Trajectory

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging, assert_logging_with_exception

# pylint: disable=protected-access

green_kubo_module = sys.modules[GreenKubo.__module__]



def _ou_velocities(n_frames, n_atoms, dt, tau, sigma, seed):
    """
    Generates an Ornstein-Uhlenbeck (AR(1)) velocity process with a
    stationary velocity auto-correlation ``sigma^2 exp(-|lag| dt/tau)``.
    """
    rng = np.random.default_rng(seed)
    phi = np.exp(-dt / tau)
    noise_scale = sigma * np.sqrt(1.0 - phi * phi)

    vel = np.zeros((n_frames, n_atoms, 3))
    vel[0] = sigma * rng.standard_normal((n_atoms, 3))
    for n in range(1, n_frames):
        vel[n] = phi * vel[n - 1] + noise_scale * rng.standard_normal(
            (n_atoms, 3)
        )

    return vel



def _make_velocity_trajectory(velocities, names=None):
    """
    Builds a trajectory from an (n_frames, n_atoms, 3) velocity array.
    """
    velocities = np.asarray(velocities, dtype=float)
    n_atoms = velocities.shape[1]

    if names is None:
        names = ["Ar"] * n_atoms

    atoms = [Atom(name) for name in names]

    systems = [
        AtomicSystem(atoms=atoms, vel=frame_velocities)
        for frame_velocities in velocities
    ]

    return Trajectory(systems)



def _write_velocity_file(path, velocities, names=None):
    """
    Writes an (n_frames, n_atoms, 3) velocity array to a .vel file.
    """
    velocities = np.asarray(velocities, dtype=float)
    n_atoms = velocities.shape[1]

    if names is None:
        names = ["Ar"] * n_atoms

    lines = []
    for frame in velocities:
        lines.append(f"{n_atoms} 100.0 100.0 100.0\n\n")
        for name, (x, y, z) in zip(names, frame):
            lines.append(f"{name} {x:.10e} {y:.10e} {z:.10e}\n")

    with open(path, "w", encoding="utf-8") as file:
        file.writelines(lines)

    return str(path)



class TestGreenKuboHelpers:

    """
    Tests for the module-level numeric helpers.
    """

    def test_constants(self):
        """
        The hardcoded unit conversion constants have the derived
        values.
        """
        assert ANGSTROM2_PER_S2_PS_TO_M2_PER_S == 1.0e-32
        assert M2_PER_S_TO_CM2_PER_S == 1.0e4

    def test_cumulative_trapezoid(self):
        """
        The cumulative trapezoid has a leading zero and reproduces the
        analytic integral of a linear integrand.
        """
        values = np.array([1.0, 2.0, 3.0, 4.0])
        result = cumulative_trapezoid(values, dx=0.5)

        assert result[0] == 0.0
        assert np.allclose(result, [0.0, 0.75, 2.0, 3.75])
        assert len(result) == len(values)

    def test_fft_matches_direct_gap_one(self):
        """
        For gap == 1 the FFT and direct estimators of the
        un-normalized velocity ACF coincide up to floating point noise.
        """
        rng = np.random.default_rng(7)
        velocities = rng.standard_normal((60, 4, 3))

        fft = velocity_acf_fft(velocities, window_size=20)
        direct = velocity_acf_direct(velocities, window_size=20, gap=1)

        assert np.allclose(fft, direct, atol=1e-10)

    def test_fft_zero_velocities_raises(self):
        """
        A trajectory with vanishing velocities has an undefined
        diffusion coefficient.
        """
        with pytest.raises(GreenKuboError, match="mean squared velocity"):
            velocity_acf_fft(np.zeros((10, 2, 3)), window_size=3)

    def test_direct_zero_velocities_raises(self):
        """
        The direct estimator raises for vanishing velocities as well.
        """
        with pytest.raises(GreenKuboError, match="mean squared velocity"):
            velocity_acf_direct(np.zeros((10, 2, 3)), window_size=3, gap=1)



class TestGreenKubo:

    """
    Tests for the GreenKubo class.
    """

    def test_defaults(self):
        """
        Without explicit parameters the class-level defaults are used
        and the window is capped to n_frames - 1.
        """
        traj = _make_velocity_trajectory(np.ones((50, 1, 3)))

        green_kubo = GreenKubo(traj, time_step=0.002)

        assert green_kubo.gap == 1
        assert green_kubo.method == "fft"
        assert green_kubo.coefficient == "diffusion"
        assert green_kubo.fit_start == 0.5
        assert green_kubo.fit_stop == 1.0
        assert green_kubo.window_size == 49
        assert green_kubo.n_blocks == 5

    def test_progress_bar_binds_config_at_call_time(self, monkeypatch):
        """
        config.with_progress_bar is set by the CLI after the module
        import, so it must be read at call time, not bound by value at
        import time.
        """
        captured = {}

        def fake_tqdm(iterable, **kwargs):
            captured.update(kwargs)
            return iterable

        monkeypatch.setattr(green_kubo_module, "tqdm", fake_tqdm)

        rng = np.random.default_rng(42)
        velocities = rng.standard_normal((20, 2, 3))

        monkeypatch.setattr(config, "with_progress_bar", False)
        GreenKubo(
            _make_velocity_trajectory(velocities),
            time_step=0.002,
            window_size=5,
        ).run()
        assert captured["disable"] is True

        captured.clear()

        monkeypatch.setattr(config, "with_progress_bar", True)
        GreenKubo(
            _make_velocity_trajectory(velocities),
            time_step=0.002,
            window_size=5,
        ).run()
        assert captured["disable"] is False

    def test_run_shapes_and_units(self):
        """
        The run returns the lag-time axis, the un-normalized velocity
        ACF and the running diffusion coefficient with consistent
        shapes and a leading zero of the running integral.
        """
        rng = np.random.default_rng(1)
        velocities = rng.standard_normal((100, 3, 3))

        green_kubo = GreenKubo(
            _make_velocity_trajectory(velocities),
            time_step=0.01,
            window_size=30,
        )
        lag_times, cvv, d_running = green_kubo.run()

        assert len(lag_times) == 31
        assert len(cvv) == 31
        assert len(d_running) == 31
        assert np.allclose(lag_times, np.arange(31) * 0.01)
        assert cvv[0] > 0.0
        assert d_running[0] == 0.0
        # cvv[0] is the aggregate mean squared velocity per atom
        expected_c0 = np.mean(np.sum(velocities ** 2, axis=2))
        assert np.isclose(cvv[0], expected_c0)

    def test_fft_matches_direct(self):
        """
        The class-level fft and direct estimators (gap == 1) give the
        same running diffusion coefficient up to floating point noise.
        """
        rng = np.random.default_rng(4711)
        velocities = rng.standard_normal((80, 4, 3)) * 1.0e12
        traj = _make_velocity_trajectory(velocities)

        fft = GreenKubo(traj, time_step=0.002, window_size=25, method="fft")
        direct = GreenKubo(
            traj,
            time_step=0.002,
            window_size=25,
            method="direct",
        )

        _, cvv_fft, d_fft = fft.run()
        _, cvv_direct, d_direct = direct.run()

        assert np.allclose(cvv_fft, cvv_direct, rtol=1e-10, atol=0.0)
        assert np.allclose(d_fft, d_direct, rtol=1e-10, atol=0.0)
        assert np.isclose(
            fft.diffusion_coefficient,
            direct.diffusion_coefficient,
            rtol=1e-10,
        )

    def test_diffusion_coefficient_units(self):
        """
        The cm^2/s properties are the m^2/s values scaled by 1e4.
        """
        rng = np.random.default_rng(2)
        velocities = rng.standard_normal((60, 2, 3)) * 1.0e12

        green_kubo = GreenKubo(
            _make_velocity_trajectory(velocities),
            time_step=0.002,
            window_size=20,
        )
        green_kubo.run()

        assert np.isclose(
            green_kubo.diffusion_coefficient_cm2_per_s,
            green_kubo.diffusion_coefficient * 1.0e4,
        )
        assert np.isclose(
            green_kubo.diffusion_coefficient_stderr_cm2_per_s,
            green_kubo.diffusion_coefficient_stderr * 1.0e4,
        )

    def test_target_selection(self):
        """
        With a target selection only the selected atoms contribute to
        the velocity ACF.
        """
        rng = np.random.default_rng(99)
        velocities = rng.standard_normal((40, 4, 3))

        selected = GreenKubo(
            _make_velocity_trajectory(velocities, names=["O", "H", "H", "O"]),
            time_step=0.002,
            window_size=10,
            target_species="O",
        )
        _, cvv_selected, _ = selected.run()

        reference = velocity_acf_fft(velocities[:, [0, 3]], window_size=10)

        assert len(selected.target_indices) == 2
        assert np.allclose(cvv_selected, reference, atol=1e-10)

    def test_raw_reader_fast_path(self, tmp_path):
        """
        A velocity TrajectoryReader uses the raw fast path and
        reproduces the in-memory result of the same (float32-rounded)
        velocities.
        """
        rng = np.random.default_rng(321)
        velocities = (rng.standard_normal((50, 3, 3)) * 1.0e12).astype(
            np.float32
        ).astype(np.float64)

        filename = _write_velocity_file(
            tmp_path / "traj.vel",
            velocities,
        )

        from_reader = GreenKubo(
            TrajectoryReader(filename),
            time_step=0.002,
            window_size=15,
        )
        from_memory = GreenKubo(
            _make_velocity_trajectory(velocities),
            time_step=0.002,
            window_size=15,
        )

        _, cvv_reader, d_reader = from_reader.run()
        _, cvv_memory, d_memory = from_memory.run()

        assert from_reader._raw_reader is not None
        assert np.allclose(cvv_reader, cvv_memory, rtol=1e-6, atol=0.0)
        assert np.allclose(d_reader, d_memory, rtol=1e-6, atol=0.0)

    def test_time_step_not_positive(self, caplog):
        """
        A non-positive time step raises a GreenKuboError.
        """
        traj = _make_velocity_trajectory(np.ones((20, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="GreenKubo",
            logging_level="ERROR",
            message_to_test="The time_step must be a positive real number.",
            exception=GreenKuboError,
            function=GreenKubo,
            traj=traj,
            time_step=0.0,
            window_size=5,
        )

    def test_unknown_method(self, caplog):
        """
        An unknown estimator method raises a GreenKuboError.
        """
        traj = _make_velocity_trajectory(np.ones((20, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="GreenKubo",
            logging_level="ERROR",
            message_to_test=(
                "Unknown method 'foo'. Possible methods are: fft, direct."
            ),
            exception=GreenKuboError,
            function=GreenKubo,
            traj=traj,
            time_step=0.002,
            window_size=5,
            method="foo",
        )

    def test_unknown_coefficient(self, caplog):
        """
        An unsupported transport coefficient raises a GreenKuboError.
        """
        traj = _make_velocity_trajectory(np.ones((20, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="GreenKubo",
            logging_level="ERROR",
            message_to_test=(
                "Unknown transport coefficient 'viscosity'. Currently "
                "supported coefficients are: diffusion."
            ),
            exception=GreenKuboError,
            function=GreenKubo,
            traj=traj,
            time_step=0.002,
            window_size=5,
            coefficient="viscosity",
        )

    def test_invalid_fit_window(self, caplog):
        """
        A fit window with fit_stop <= fit_start raises a GreenKuboError.
        """
        traj = _make_velocity_trajectory(np.ones((20, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="GreenKubo",
            logging_level="ERROR",
            message_to_test=(
                "The fit window fractions fit_start=0.8 and fit_stop=0.5 "
                "must satisfy 0 <= fit_start < fit_stop <= 1."
            ),
            exception=GreenKuboError,
            function=GreenKubo,
            traj=traj,
            time_step=0.002,
            window_size=5,
            fit_start=0.8,
            fit_stop=0.5,
        )

    def test_window_too_long(self, caplog):
        """
        A window that is not shorter than the trajectory raises a
        GreenKuboError.
        """
        traj = _make_velocity_trajectory(np.ones((20, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="GreenKubo",
            logging_level="ERROR",
            message_to_test=(
                "The trajectory contains only 20 frame(s), but at least "
                "window_size + 1 = 21 frames are needed so that every "
                "lag has at least one time origin."
            ),
            exception=GreenKuboError,
            function=GreenKubo,
            traj=traj,
            time_step=0.002,
            window_size=20,
        )

    def test_empty_trajectory(self, caplog):
        """
        An empty trajectory raises a GreenKuboError.
        """
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="GreenKubo",
            logging_level="ERROR",
            message_to_test="Trajectory cannot be of length 0.",
            exception=GreenKuboError,
            function=GreenKubo,
            traj=Trajectory(),
            time_step=0.002,
            window_size=5,
        )

    def test_empty_target_selection(self, caplog):
        """
        A target selection that selects no atoms raises a
        GreenKuboError.
        """
        traj = _make_velocity_trajectory(
            np.ones((20, 2, 3)),
            names=["O", "H"],
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="GreenKubo",
            logging_level="ERROR",
            message_to_test="The target selection does not select any atoms.",
            exception=GreenKuboError,
            function=GreenKubo,
            traj=traj,
            time_step=0.002,
            window_size=5,
            target_species="C",
        )

    def test_missing_velocities(self, caplog):
        """
        A trajectory without velocities raises a GreenKuboError at run
        time.
        """
        atoms = [Atom("Ar")]
        systems = [
            AtomicSystem(atoms=atoms, pos=np.zeros((1, 3)))
            for _ in range(20)
        ]
        green_kubo = GreenKubo(
            Trajectory(systems),
            time_step=0.002,
            window_size=5,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="GreenKubo",
            logging_level="ERROR",
            message_to_test=(
                "A frame of the velocity trajectory does not provide "
                "velocities for all 1 atoms. Please provide a velocity "
                "trajectory (e.g. .vel files)."
            ),
            exception=GreenKuboError,
            function=green_kubo.run,
        )



class TestGreenKuboBlockAveraging:

    """
    Tests for the block-averaged statistical uncertainty of the
    diffusion coefficient and the n_blocks handling.
    """

    def test_block_sem_and_plateau_spread_are_distinct(self):
        """
        The reported uncertainty is the block-averaged standard error
        (over n_blocks blocks) and is stored separately from the
        (honestly renamed) plateau spread of the running integral.
        """
        velocities = _ou_velocities(
            n_frames=2000,
            n_atoms=8,
            dt=0.002,
            tau=0.1,
            sigma=2.0e12,
            seed=4000,
        )

        green_kubo = GreenKubo(
            _make_velocity_trajectory(velocities),
            time_step=0.002,
            window_size=200,
            n_blocks=5,
        )
        green_kubo.run()

        # five per-block diffusion coefficients feed the standard error
        assert green_kubo.n_blocks == 5
        assert len(green_kubo.block_diffusion_coefficients) == 5

        expected_stderr = float(
            np.std(green_kubo.block_diffusion_coefficients, ddof=1)
            / np.sqrt(5)
        )
        assert green_kubo.diffusion_coefficient_stderr == pytest.approx(
            expected_stderr,
            rel=1e-12,
        )

        # the plateau spread is a different quantity and here clearly
        # smaller than the (correct) block standard error
        assert green_kubo.diffusion_coefficient_plateau_spread > 0.0
        assert (
            green_kubo.diffusion_coefficient_plateau_spread
            < green_kubo.diffusion_coefficient_stderr
        )

    def test_stderr_cm2_reflects_block_sem(self):
        """
        The cm^2/s uncertainty property scales the block standard error
        by 1e4.
        """
        velocities = _ou_velocities(
            n_frames=1200,
            n_atoms=6,
            dt=0.002,
            tau=0.1,
            sigma=2.0e12,
            seed=17,
        )

        green_kubo = GreenKubo(
            _make_velocity_trajectory(velocities),
            time_step=0.002,
            window_size=150,
            n_blocks=4,
        )
        green_kubo.run()

        assert np.isclose(
            green_kubo.diffusion_coefficient_stderr_cm2_per_s,
            green_kubo.diffusion_coefficient_stderr * 1.0e4,
        )

    def test_explicit_n_blocks(self):
        """
        An explicit n_blocks that fits the trajectory is used verbatim.
        """
        velocities = _ou_velocities(
            n_frames=1500,
            n_atoms=4,
            dt=0.002,
            tau=0.1,
            sigma=2.0e12,
            seed=3,
        )

        green_kubo = GreenKubo(
            _make_velocity_trajectory(velocities),
            time_step=0.002,
            window_size=100,
            n_blocks=6,
        )
        green_kubo.run()

        assert len(green_kubo.block_diffusion_coefficients) == 6

    def test_n_blocks_too_small(self, caplog):
        """
        An n_blocks smaller than 2 has no defined block standard error
        and raises a GreenKuboError.
        """
        traj = _make_velocity_trajectory(np.ones((20, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="GreenKubo",
            logging_level="ERROR",
            message_to_test=(
                "The number of blocks n_blocks=1 must be at least 2 so "
                "that a block-averaged standard error of the diffusion "
                "coefficient can be computed."
            ),
            exception=GreenKuboError,
            function=GreenKubo,
            traj=traj,
            time_step=0.002,
            window_size=5,
            n_blocks=1,
        )

    def test_too_few_frames_per_block_raises(self, caplog):
        """
        A trajectory that is long enough for the window but too short to
        form even two blocks that each cover the window raises a
        GreenKuboError at run time.
        """
        rng = np.random.default_rng(5)
        velocities = rng.standard_normal((12, 2, 3)) * 1.0e12

        green_kubo = GreenKubo(
            _make_velocity_trajectory(velocities),
            time_step=0.002,
            window_size=8,
            n_blocks=5,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="GreenKubo",
            logging_level="ERROR",
            message_to_test=(
                "The trajectory (12 frames, correlation window 8) is too "
                "short to be split into 5 blocks that each cover the "
                "correlation window. At least 2 blocks longer than "
                "window_size = 8 frames are needed, but only 1 fit(s). "
                "Reduce n_blocks or the window size, or provide a longer "
                "trajectory."
            ),
            exception=GreenKuboError,
            function=green_kubo.run,
        )

    def test_n_blocks_clamped_with_warning(self, caplog):
        """
        A requested n_blocks that is too large for the window is clamped
        down to the largest feasible number of blocks and a warning is
        logged.
        """
        rng = np.random.default_rng(6)
        velocities = rng.standard_normal((100, 3, 3)) * 1.0e12

        green_kubo = GreenKubo(
            _make_velocity_trajectory(velocities),
            time_step=0.002,
            window_size=30,
            n_blocks=5,
        )

        assert_logging(
            caplog=caplog,
            logging_name="GreenKubo",
            logging_level="WARNING",
            message_to_test=(
                "The requested number of blocks n_blocks=5 is too large "
                "for a trajectory of 100 frames with correlation window "
                "30; each block must be longer than the window. Clamping "
                "to 3 blocks."
            ),
            function=green_kubo.run,
        )

        # the block average actually used the clamped number of blocks
        assert len(green_kubo.block_diffusion_coefficients) == 3
