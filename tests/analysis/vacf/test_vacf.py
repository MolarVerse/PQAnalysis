"""
Tests for the VACF class.
"""

import sys

import numpy as np
import pytest

from PQAnalysis import config
from PQAnalysis.analysis.vacf import VACF
from PQAnalysis.analysis.vacf.exceptions import VACFError
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.traj import Trajectory

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception

# pylint: disable=protected-access

vacf_module = sys.modules[VACF.__module__]



def _make_velocity_trajectory(velocities, names=None):
    """
    Builds a trajectory from an (n_frames, n_atoms, 3) velocity array.
    """
    velocities = np.asarray(velocities, dtype=float)
    n_atoms = velocities.shape[1]

    if names is None:
        names = ["O"] * n_atoms

    atoms = [Atom(name) for name in names]

    systems = [
        AtomicSystem(atoms=atoms, vel=frame_velocities)
        for frame_velocities in velocities
    ]

    return Trajectory(systems)



def _make_charge_trajectory(charges, names=None):
    """
    Builds a charge trajectory from an (n_frames, n_atoms) charge array.
    """
    charges = np.asarray(charges, dtype=float)
    n_atoms = charges.shape[1]

    if names is None:
        names = ["O"] * n_atoms

    atoms = [Atom(name) for name in names]

    systems = [
        AtomicSystem(atoms=atoms, charges=frame_charges)
        for frame_charges in charges
    ]

    return Trajectory(systems)



def _reference_vacf(velocities, window_size, gap):
    """
    Brute-force emulation of the legacy FreqCalc estimator.
    """
    velocities = np.asarray(velocities, dtype=float)
    n_frames = velocities.shape[0]

    stop_frame = (n_frames - window_size) // gap * gap

    corr = np.zeros(window_size + 1)
    n_origins = 0

    for origin in range(gap, stop_frame + 1, gap):
        origin_velocities = velocities[origin - 1]
        norm = np.sum(origin_velocities * origin_velocities)
        n_origins += 1

        for lag in range(window_size + 1):
            corr[lag] += np.sum(
                velocities[origin - 1 + lag] * origin_velocities
            ) / norm

    return corr / n_origins



def _reference_fft_vacf(velocities, window_size):
    """
    Brute-force emulation of the Wiener-Khinchin (fft) estimator.
    """
    velocities = np.asarray(velocities, dtype=float)
    n_frames = velocities.shape[0]

    raw = np.zeros(window_size + 1)

    for lag in range(window_size + 1):
        scalar = 0.0
        for frame in range(n_frames - lag):
            scalar += np.sum(
                velocities[frame] * velocities[frame + lag]
            )
        raw[lag] = scalar / (n_frames - lag)

    return raw / raw[0]



#: Velocities of the n_frames == window_size boundary-case parity
#: check (20 frames, 2 atoms). The corresponding QMCFC velocity file
#: (values printed with 7 decimals) was run through the recompiled
#: legacy FreqCalc binary with window = 20, gap = 1, time_step = 0.1
#: and target_atoms = 1-2.
BOUNDARY_VELOCITIES = np.array(
    [
        [[-0.7931225, 0.2405713, -1.8963263],
            [1.3957717, 0.6382947, -0.2920475]],
        [[-0.3119493, 0.3038354, -0.2676603],
            [-0.2259089, 0.7200678, 0.5147052]],
        [[-0.0641279, -0.0854766, 0.1609163],
            [-0.6140184, -0.4037503, 0.5482602]],
        [[-0.1304828, -1.3744262, -0.4772787],
            [0.6566216, -0.2322828, -0.1487328]],
        [[0.6418366, 1.8246103, -0.7131887],
            [1.3482068, -1.2300128, 0.1749776]],
        [[-1.1695295, 1.3514582, 0.8339230],
            [1.1377168, -0.8855332, 0.6845550]],
        [[-0.5190133, -0.4573853, 0.5065371],
            [0.8767182, 0.2044202, -0.6279871]],
        [[-0.8258165, 1.4443170, 0.5939459],
            [0.7197279, 2.1834881, -0.8158629]],
        [[2.5594578, 3.1509081, 1.6184131],
            [0.8271065, -0.6638220, 0.9944856]],
        [[-0.4426921, -0.0216508, -0.2904332],
            [0.2838300, 1.2880777, -0.5555841]],
        [[-0.9854038, -1.0029566, -0.9682773],
            [-1.4311030, -0.9129341, 1.2932656]],
        [[-0.5933120, 0.2571077, -1.2168940],
            [0.1696506, -1.7407814, -0.6987559]],
        [[2.2545389, -0.5829697, 1.1199764],
            [0.4550541, -0.1529858, -0.6521068]],
        [[1.2869094, -0.1773835, 1.5274245],
            [-0.7190790, 0.0573390, 0.4654991]],
        [[0.3731602, -1.2338015, -0.6640627],
            [-0.1959802, -0.8536993, 0.6773251]],
        [[0.5880192, -1.9570845, -1.8052541],
            [-1.2815633, 0.1172595, 2.0331725]],
        [[-0.3823565, 0.2506505, -1.0631117],
            [-1.0468362, -1.9572278, -0.0283442]],
        [[0.9472167, -0.3567143, 1.3964531],
            [0.1978549, -0.0364125, 0.5172539]],
        [[0.4873776, 1.1478064, -0.8019543],
            [-2.2880492, 0.1147467, -0.6116752]],
        [[-0.0271772, 1.6643000, -1.1028418],
            [0.7647973, 0.9458675, 0.4607362]],
    ]
)

#: The VACF of BOUNDARY_VELOCITIES as printed (8 decimals) by the
#: recompiled legacy FreqCalc binary: a single origin is spawned at
#: frame 1 (stop_frame reset 0 -> 1) and the final lag bin stays zero.
BOUNDARY_LEGACY_VACF = np.array(
    [
        1.00000000, 0.12225998, -0.23047930,
        0.22153690, 0.34621604, 0.07349218,
        0.13067321, 0.37369523, -0.58012556,
        0.33868165, -0.08647947, 0.32269875,
        -0.49439806, -0.75296970, -0.05201860,
        0.02654049, -0.04797735, -0.50318541,
        -0.22775122, 0.60232549, 0.00000000,
    ]
)



class TestVACF:

    """
    Tests for the VACF class.
    """

    def test_window_size_default(self):
        """
        Without an explicit window_size the class-level default of
        1000 frames is used.
        """
        traj = _make_velocity_trajectory(np.ones((1000, 1, 3)))

        vacf = VACF(traj, time_step=0.1)

        assert vacf.window_size == 1000
        assert vacf.window_size == VACF._window_size_default

    def test_progress_bar_binds_config_at_call_time(self, monkeypatch):
        """
        config.with_progress_bar is set by the CLI after the module
        import, so it must be read at call time, not bound by value
        at import time.
        """
        captured = {}

        def fake_tqdm(iterable, **kwargs):
            captured.update(kwargs)
            return iterable

        monkeypatch.setattr(vacf_module, "tqdm", fake_tqdm)

        rng = np.random.default_rng(42)
        velocities = rng.standard_normal((6, 2, 3))

        monkeypatch.setattr(config, "with_progress_bar", False)
        VACF(
            _make_velocity_trajectory(velocities),
            window_size=2,
            time_step=0.1,
        ).run()
        assert captured["disable"] is True

        captured.clear()

        monkeypatch.setattr(config, "with_progress_bar", True)
        VACF(
            _make_velocity_trajectory(velocities),
            window_size=2,
            time_step=0.1,
        ).run()
        assert captured["disable"] is False

    def test_c0_is_one(self):
        """
        The normalized VACF starts at exactly one.
        """
        rng = np.random.default_rng(42)
        traj = _make_velocity_trajectory(
            rng.standard_normal((50, 3, 3))
        )

        vacf = VACF(traj, window_size=10, time_step=0.5, gap=2)
        time, correlation = vacf.run()

        assert correlation[0] == 1.0
        assert np.allclose(time, np.arange(11) * 0.5, rtol=1e-14)
        assert len(correlation) == 11

    def test_matches_brute_force_reference(self):
        """
        The vectorized sliding-origin estimator matches a brute-force
        emulation of the legacy FreqCalc algorithm, including the
        drain phase and the origin counting.
        """
        rng = np.random.default_rng(4711)
        velocities = rng.standard_normal((43, 4, 3))
        traj = _make_velocity_trajectory(velocities)

        vacf = VACF(traj, window_size=12, time_step=0.1, gap=3)
        _, correlation = vacf.run()

        reference = _reference_vacf(velocities, window_size=12, gap=3)

        assert vacf.n_origins == (43 - 12) // 3
        assert np.allclose(correlation, reference, atol=1e-12)

    def test_matches_brute_force_reference_gap_one(self):
        """
        The estimator matches the brute-force reference for gap 1.
        """
        rng = np.random.default_rng(1234)
        velocities = rng.standard_normal((30, 2, 3))
        traj = _make_velocity_trajectory(velocities)

        vacf = VACF(traj, window_size=8, time_step=0.1)
        _, correlation = vacf.run()

        reference = _reference_vacf(velocities, window_size=8, gap=1)

        assert np.allclose(correlation, reference, atol=1e-12)

    def test_window_equals_n_frames_legacy_parity(self):
        """
        The boundary case n_frames == window_size with gap == 1 spawns
        a single origin at frame 1 (legacy FreqCalc stop_frame reset)
        and matches the recompiled legacy binary; only the final lag
        bin stays zero.
        """
        traj = _make_velocity_trajectory(
            BOUNDARY_VELOCITIES,
            names=["O", "H"],
        )

        vacf = VACF(traj, window_size=20, time_step=0.1, gap=1)
        _, correlation = vacf.run()

        assert vacf.stop_frame == 1
        assert vacf.n_origins == 1
        assert correlation[0] == 1.0
        assert correlation[-1] == 0.0
        # tolerance: the legacy binary prints with 8 decimals
        assert np.allclose(
            correlation,
            BOUNDARY_LEGACY_VACF,
            atol=5e-8,
        )

    def test_window_equals_n_frames_gap_two_rejected(self, caplog):
        """
        The boundary case n_frames == window_size is only valid for
        gap == 1; for larger gaps no origin can spawn (the legacy tool
        would divide by zero) and a VACFError is raised.
        """
        traj = _make_velocity_trajectory(np.ones((20, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "The trajectory contains only 20 frame(s), but at "
                "least window_size + gap = 22 frames are needed to "
                "place a single time origin (or exactly window_size "
                "frames for the direct method with gap == 1)."
            ),
            exception=VACFError,
            function=VACF,
            traj=traj,
            window_size=20,
            time_step=0.1,
            gap=2,
        )

    def test_window_equals_n_frames_fft_rejected(self, caplog):
        """
        The fft estimator needs at least one origin for every lag and
        therefore still rejects n_frames == window_size with gap == 1.
        """
        traj = _make_velocity_trajectory(np.ones((20, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "The trajectory contains only 20 frame(s), but at "
                "least window_size + gap = 21 frames are needed to "
                "place a single time origin (or exactly window_size "
                "frames for the direct method with gap == 1)."
            ),
            exception=VACFError,
            function=VACF,
            traj=traj,
            window_size=20,
            time_step=0.1,
            method="fft",
        )

    @pytest.mark.parametrize("method", ["direct", "fft"])
    def test_single_cosine_analytic(self, method):
        """
        For velocities v(t) = cos(2 pi nu t + phi) with evenly
        distributed phases the VACF is the analytic cosine
        cos(2 pi nu lag dt) up to floating point noise.

        This is an implementation-independent physical check for both
        estimators. (The direct estimator with gap = 1 and the fft
        estimator do not coincide on general data - they use different
        origin sets and normalizations - so the analytic case is the
        common ground truth.)
        """
        n_atoms = 8
        n_frames = 200
        time_step = 0.002
        frequency = 25.0  # ps^-1

        time = np.arange(n_frames) * time_step
        phases = (
            2.0 * np.pi * np.arange(3 * n_atoms) / (3 * n_atoms)
        ).reshape(n_atoms, 3)

        velocities = np.cos(
            2.0 * np.pi * frequency * time[:, None, None] + phases[None]
        )

        traj = _make_velocity_trajectory(velocities)
        vacf = VACF(
            traj,
            window_size=100,
            time_step=time_step,
            gap=5,
            method=method,
        )
        lag_time, correlation = vacf.run()

        analytic = np.cos(2.0 * np.pi * frequency * lag_time)

        assert np.allclose(correlation, analytic, atol=1e-10)

    def test_fft_matches_brute_force_reference(self):
        """
        The fft method matches a brute-force emulation of the
        denser-origin Wiener-Khinchin estimator.
        """
        rng = np.random.default_rng(7)
        velocities = rng.standard_normal((32, 3, 3))
        traj = _make_velocity_trajectory(velocities)

        vacf = VACF(traj, window_size=10, time_step=0.1, method="fft")
        _, correlation = vacf.run()

        reference = _reference_fft_vacf(velocities, window_size=10)

        assert correlation[0] == 1.0
        assert vacf.n_origins == 32
        assert np.allclose(correlation, reference, atol=1e-12)

    def test_static_charges_weighting(self):
        """
        The charge-flux mode with static charges is identical to the
        plain VACF of the charge-scaled velocities.
        """
        rng = np.random.default_rng(2024)
        velocities = rng.standard_normal((25, 3, 3))
        charges = np.array([-0.8, 0.4, 0.4])

        flux = VACF(
            _make_velocity_trajectory(velocities),
            window_size=6,
            time_step=0.1,
            gap=2,
            charges=charges,
        )
        _, flux_correlation = flux.run()

        scaled = VACF(
            _make_velocity_trajectory(velocities * charges[None, :, None]),
            window_size=6,
            time_step=0.1,
            gap=2,
        )
        _, scaled_correlation = scaled.run()

        assert flux.flux
        assert not scaled.flux
        assert np.allclose(flux_correlation, scaled_correlation, atol=1e-14)

    def test_charge_trajectory_weighting(self):
        """
        The charge-flux mode with a charge trajectory is identical to
        the plain VACF of the frame-wise charge-scaled velocities.
        """
        rng = np.random.default_rng(2025)
        velocities = rng.standard_normal((25, 3, 3))
        charges = 0.5 + 0.1 * rng.standard_normal((25, 3))

        flux = VACF(
            _make_velocity_trajectory(velocities),
            window_size=6,
            time_step=0.1,
            gap=2,
            charge_traj=_make_charge_trajectory(charges),
        )
        _, flux_correlation = flux.run()

        scaled = VACF(
            _make_velocity_trajectory(velocities * charges[:, :, None]),
            window_size=6,
            time_step=0.1,
            gap=2,
        )
        _, scaled_correlation = scaled.run()

        assert flux.flux
        assert np.allclose(flux_correlation, scaled_correlation, atol=1e-14)

    def test_target_selection(self):
        """
        With a target selection only the selected atoms contribute.
        """
        rng = np.random.default_rng(99)
        velocities = rng.standard_normal((25, 4, 3))

        selected = VACF(
            _make_velocity_trajectory(
                velocities,
                names=["O", "H", "H", "O"],
            ),
            window_size=6,
            time_step=0.1,
            gap=2,
            target_species="O",
        )
        _, selected_correlation = selected.run()

        reference = _reference_vacf(
            velocities[:, [0, 3]],
            window_size=6,
            gap=2,
        )

        assert np.allclose(selected_correlation, reference, atol=1e-12)

    def test_unknown_method(self, caplog):
        """
        An unknown estimator method raises a VACFError.
        """
        traj = _make_velocity_trajectory(np.ones((20, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "Unknown method 'foo'. Possible methods are: direct, fft."
            ),
            exception=VACFError,
            function=VACF,
            traj=traj,
            window_size=5,
            time_step=0.1,
            method="foo",
        )

    def test_both_charge_sources(self, caplog):
        """
        Static charges and a charge trajectory cannot be combined.
        """
        traj = _make_velocity_trajectory(np.ones((20, 1, 3)))
        charge_traj = _make_charge_trajectory(np.ones((20, 1)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "Only one charge source can be used for the charge-flux "
                "mode: either static charges or a charge trajectory."
            ),
            exception=VACFError,
            function=VACF,
            traj=traj,
            window_size=5,
            time_step=0.1,
            charges=np.ones(1),
            charge_traj=charge_traj,
        )

    def test_window_size_not_multiple_of_gap(self, caplog):
        """
        The window size has to be an integer multiple of the gap.
        """
        traj = _make_velocity_trajectory(np.ones((30, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "The window_size 10 must be an integer multiple of the "
                "gap 3 for the sliding-origin machinery."
            ),
            exception=VACFError,
            function=VACF,
            traj=traj,
            window_size=10,
            time_step=0.1,
            gap=3,
        )

    def test_trajectory_too_short(self, caplog):
        """
        The trajectory has to accommodate at least one time origin.
        """
        traj = _make_velocity_trajectory(np.ones((10, 1, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "The trajectory contains only 10 frame(s), but at least "
                "window_size + gap = 12 frames are needed to place a "
                "single time origin (or exactly window_size frames for "
                "the direct method with gap == 1)."
            ),
            exception=VACFError,
            function=VACF,
            traj=traj,
            window_size=10,
            time_step=0.1,
            gap=2,
        )

    def test_empty_target_selection(self, caplog):
        """
        A target selection that does not select any atoms raises a
        VACFError instead of a misleading zero-norm error at run time.
        """
        traj = _make_velocity_trajectory(
            np.ones((20, 2, 3)),
            names=["O", "H"],
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "The target selection does not select any atoms."
            ),
            exception=VACFError,
            function=VACF,
            traj=traj,
            window_size=5,
            time_step=0.1,
            target_species="C",
        )

    def test_empty_trajectory(self, caplog):
        """
        An empty trajectory raises a VACFError.
        """
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test="Trajectory cannot be of length 0.",
            exception=VACFError,
            function=VACF,
            traj=Trajectory(),
            window_size=5,
            time_step=0.1,
        )

    def test_wrong_number_of_static_charges(self, caplog):
        """
        The number of static charges has to match the number of atoms.
        """
        traj = _make_velocity_trajectory(np.ones((20, 3, 3)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "The number of static charges 2 does not match the "
                "number of atoms 3 of the system."
            ),
            exception=VACFError,
            function=VACF,
            traj=traj,
            window_size=5,
            time_step=0.1,
            charges=np.ones(2),
        )

    def test_charge_trajectory_not_in_lockstep(self, caplog):
        """
        The charge trajectory has to provide exactly one frame per
        velocity frame.
        """
        traj = _make_velocity_trajectory(np.ones((20, 1, 3)))
        charge_traj = _make_charge_trajectory(np.ones((19, 1)))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "The charge trajectory contains 19 frame(s), but the "
                "velocity trajectory contains 20 frame(s). Both "
                "trajectories have to be in lockstep."
            ),
            exception=VACFError,
            function=VACF,
            traj=traj,
            window_size=5,
            time_step=0.1,
            charge_traj=charge_traj,
        )

    def test_zero_norm_origin(self, caplog):
        """
        A time origin with a vanishing aggregate squared velocity norm
        raises a VACFError.
        """
        traj = _make_velocity_trajectory(np.zeros((20, 1, 3)))
        vacf = VACF(traj, window_size=5, time_step=0.1)

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "The aggregate squared velocity norm of the time origin "
                "at frame 1 is zero. The normalized VACF is not defined."
            ),
            exception=VACFError,
            function=vacf.run,
        )

    def test_missing_velocities(self, caplog):
        """
        A trajectory without velocities raises a VACFError at run time.
        """
        atoms = [Atom("O")]
        systems = [
            AtomicSystem(atoms=atoms, pos=np.zeros((1, 3)))
            for _ in range(20)
        ]
        vacf = VACF(Trajectory(systems), window_size=5, time_step=0.1)

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "A frame of the velocity trajectory does not provide "
                "velocities for all 1 atoms. Please provide a velocity "
                "trajectory (e.g. .vel files)."
            ),
            exception=VACFError,
            function=vacf.run,
        )
