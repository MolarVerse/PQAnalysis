"""
Physical validation of the Green-Kubo self-diffusion coefficient.

The tests cover four independent checks:

(a) ANALYTIC:  a synthetic Ornstein-Uhlenbeck (AR(1)) velocity process
    whose velocity auto-correlation is ``sigma^2 exp(-t/tau)`` must
    yield the known plateau ``D = sigma^2 * tau`` (in the module's unit
    system).
(b) NUMPY REFERENCE:  an independent numpy implementation of the whole
    pipeline (FFT ACF + cumulative trapezoid + prefactor) must match
    the module to ~1e-10 on the same synthetic data.
(c) MSD CROSS-CHECK:  on a trajectory that has both positions (.xyz)
    and consistent velocities (.vel), the Einstein-MSD diffusion
    coefficient and the Green-Kubo velocity-ACF diffusion coefficient
    must agree within statistics, validating the whole units chain.
(d) INVARIANCE:  fft vs direct estimator, stability against the window
    size and against sub-sampling dt -> 2 dt.
"""

import numpy as np
import pytest

from PQAnalysis import config
from PQAnalysis.analysis.green_kubo import (
    GreenKubo,
    ANGSTROM2_PER_S2_PS_TO_M2_PER_S,
)
from PQAnalysis.analysis.msd import MSD
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.traj import Trajectory

from .. import pytestmark  # pylint: disable=unused-import

# pylint: disable=protected-access

#: The analytic diffusion coefficient of the generated cross-check
#: trajectory: D = sigma^2 * tau = (2 Angstrom/ps)^2 * 0.1 ps
#: = 0.4 Angstrom^2/ps = 4e-9 m^2/s.
D_ANALYTIC_M2_PER_S = 4.0e-9



def _ou_velocities(n_frames, n_atoms, dt, tau, sigma, seed):
    """
    Generates an Ornstein-Uhlenbeck (AR(1)) velocity process with a
    per-component stationary standard deviation ``sigma`` and a
    correlation time ``tau``. The stationary velocity auto-correlation
    is ``sigma^2 exp(-|lag| dt / tau)``.
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



def _make_velocity_trajectory(velocities):
    """
    Builds a trajectory from an (n_frames, n_atoms, 3) velocity array.
    """
    atoms = [Atom("Ar") for _ in range(velocities.shape[1])]

    return Trajectory([
        AtomicSystem(atoms=atoms, vel=frame)
        for frame in velocities
    ])



def _numpy_reference(velocities, window_size, dt, fit_start, fit_stop):
    """
    Independent numpy implementation of the whole Green-Kubo pipeline.
    """
    n_frames, n_sel, _ = velocities.shape

    n_fft = 2 * n_frames
    spectrum = np.fft.rfft(velocities, n=n_fft, axis=0)
    autocorr = np.fft.irfft(
        spectrum * np.conj(spectrum),
        n=n_fft,
        axis=0,
    )[:window_size + 1].real

    raw = autocorr.sum(axis=(1, 2))
    counts = n_frames - np.arange(window_size + 1)
    cvv = raw / counts / n_sel

    increments = 0.5 * (cvv[1:] + cvv[:-1]) * dt
    integral = np.concatenate(([0.0], np.cumsum(increments)))
    d_running = ANGSTROM2_PER_S2_PS_TO_M2_PER_S / 3.0 * integral

    start = int(round(fit_start * window_size))
    stop = max(int(round(fit_stop * window_size)), start + 1)
    plateau = d_running[start:stop + 1]

    return cvv, d_running, float(np.mean(plateau))



class TestGreenKuboAnalytic:

    """
    (a) analytic plateau and (b) numpy reference checks.
    """

    def test_analytic_plateau(self):
        """
        The Green-Kubo plateau reproduces the analytic OU diffusion
        coefficient D = sigma^2 * tau within a few percent.
        """
        dt = 0.002
        tau = 0.1
        sigma = 2.0e12  # Angstrom / s

        velocities = _ou_velocities(
            n_frames=4000,
            n_atoms=400,
            dt=dt,
            tau=tau,
            sigma=sigma,
            seed=2024,
        )

        green_kubo = GreenKubo(
            _make_velocity_trajectory(velocities),
            time_step=dt,
            window_size=400,
        )
        green_kubo.run()

        # D = sigma^2 * tau in (Angstrom/s)^2 * ps, converted to m^2/s
        d_analytic = sigma ** 2 * tau * ANGSTROM2_PER_S2_PS_TO_M2_PER_S

        assert green_kubo.diffusion_coefficient == pytest.approx(
            d_analytic,
            rel=0.05,
        )

    def test_numpy_reference(self):
        """
        The module reproduces an independent numpy implementation of the
        full pipeline to ~1e-10 (relative) on the same float64 data.
        """
        dt = 0.002
        velocities = _ou_velocities(
            n_frames=1500,
            n_atoms=16,
            dt=dt,
            tau=0.1,
            sigma=2.0e12,
            seed=7,
        )

        green_kubo = GreenKubo(
            _make_velocity_trajectory(velocities),
            time_step=dt,
            window_size=400,
            fit_start=0.5,
            fit_stop=1.0,
        )
        _, cvv, d_running = green_kubo.run()

        ref_cvv, ref_running, ref_d = _numpy_reference(
            velocities,
            window_size=400,
            dt=dt,
            fit_start=0.5,
            fit_stop=1.0,
        )

        assert np.allclose(cvv, ref_cvv, rtol=1e-10, atol=0.0)
        assert np.allclose(d_running, ref_running, rtol=1e-10, atol=0.0)
        assert green_kubo.diffusion_coefficient == pytest.approx(
            ref_d,
            rel=1e-10,
        )



class TestGreenKuboMSDCrossCheck:

    """
    (c) the key physical cross-check between the Einstein-MSD and the
    Green-Kubo self-diffusion coefficient on a consistent xyz+vel
    trajectory.
    """

    @pytest.mark.parametrize("example_dir", ["green_kubo"], indirect=False)
    def test_msd_cross_check(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        The Green-Kubo velocity-ACF diffusion coefficient agrees with
        the Einstein-MSD diffusion coefficient of the same dynamics
        within ~5 %, and both agree with the analytic value within
        ~8 %. The tolerance covers the finite-trajectory statistical
        error and the different finite-time biases of the two
        estimators (the Einstein fit uses the MSD tail slope, the
        Green-Kubo estimate the velocity-ACF plateau); in expectation
        both estimators are identical for the underlying discrete
        process x_{n+1} = x_n + v_n dt.
        """
        green_kubo = GreenKubo(
            TrajectoryReader("traj.vel"),
            time_step=0.002,
            window_size=800,
        )
        green_kubo.run()

        msd = MSD(
            TrajectoryReader("traj.xyz"),
            target_species="Ar",
            window=800,
            gap=10,
            time_step=0.002,
            fit_window=160,
        )
        msd.run()

        d_gk = green_kubo.diffusion_coefficient
        d_einstein = msd.fit_results["total"].diffusion_coefficient

        assert d_gk == pytest.approx(d_einstein, rel=0.05)
        assert d_gk == pytest.approx(D_ANALYTIC_M2_PER_S, rel=0.08)
        assert d_einstein == pytest.approx(D_ANALYTIC_M2_PER_S, rel=0.08)



class TestGreenKuboBlockError:

    """
    (e) the block-averaged standard error is a faithful estimate of the
    true run-to-run scatter of the diffusion coefficient, in contrast to
    the (old) plateau spread which severely understates it.

    The check generates many independent OU trajectories, measures the
    empirical scatter of the full-trajectory diffusion coefficient
    across them (the ground-truth statistical error) and compares it to
    the two per-trajectory uncertainty estimates. The test is fully
    deterministic (fixed seeds) and uses generous factor-of-three
    tolerances so that it is robust and not flaky.
    """

    def test_block_sem_matches_true_scatter(self, monkeypatch):
        """
        Averaged over independent trajectories the block standard error
        estimates the true run-to-run scatter of D to within a factor
        of ~3 (measured ratio ~0.9), whereas the plateau spread of the
        running integral is far smaller than the true scatter (measured
        ratio ~0.3, i.e. it understates the error by ~3x).
        """
        monkeypatch.setattr(config, "with_progress_bar", False)

        n_trajectories = 24
        dt = 0.002
        tau = 0.1
        sigma = 2.0e12
        window_size = 200
        n_blocks = 5

        diffusion_coefficients = []
        block_sems = []
        plateau_spreads = []

        for index in range(n_trajectories):
            velocities = _ou_velocities(
                n_frames=2000,
                n_atoms=8,
                dt=dt,
                tau=tau,
                sigma=sigma,
                seed=4000 + index,
            )

            green_kubo = GreenKubo(
                _make_velocity_trajectory(velocities),
                time_step=dt,
                window_size=window_size,
                n_blocks=n_blocks,
            )
            green_kubo.run()

            diffusion_coefficients.append(green_kubo.diffusion_coefficient)
            block_sems.append(green_kubo.diffusion_coefficient_stderr)
            plateau_spreads.append(
                green_kubo.diffusion_coefficient_plateau_spread
            )

        # ground truth: the run-to-run scatter of the full-trajectory D
        true_scatter = float(np.std(diffusion_coefficients, ddof=1))
        mean_block_sem = float(np.mean(block_sems))
        mean_plateau_spread = float(np.mean(plateau_spreads))

        # the NEW block standard error is a reasonable estimate of the
        # true scatter (within a factor of ~3, actually ~0.9)
        assert 0.3 * true_scatter <= mean_block_sem <= 3.0 * true_scatter

        # the OLD plateau spread severely understates the true scatter
        # (here by ~3x); this is exactly the reported defect
        assert mean_plateau_spread < 0.6 * true_scatter

        # and the block standard error is much closer to the truth than
        # the plateau spread ever was
        assert abs(mean_block_sem - true_scatter) < abs(
            mean_plateau_spread - true_scatter
        )



class TestGreenKuboInvariance:

    """
    (d) estimator, window-size and sub-sampling invariance.
    """

    def test_fft_matches_direct(self):
        """
        The fft and direct estimators (gap == 1) give the same plateau
        diffusion coefficient up to floating point noise.
        """
        velocities = _ou_velocities(
            n_frames=2000,
            n_atoms=32,
            dt=0.002,
            tau=0.1,
            sigma=2.0e12,
            seed=11,
        )
        traj = _make_velocity_trajectory(velocities)

        fft = GreenKubo(traj, time_step=0.002, window_size=500, method="fft")
        direct = GreenKubo(
            traj,
            time_step=0.002,
            window_size=500,
            method="direct",
        )
        fft.run()
        direct.run()

        assert fft.diffusion_coefficient == pytest.approx(
            direct.diffusion_coefficient,
            rel=1e-9,
        )

    def test_window_size_stability(self):
        """
        Once the correlation window covers the velocity-ACF decay, the
        plateau diffusion coefficient is stable against a further
        increase of the window size.
        """
        velocities = _ou_velocities(
            n_frames=4000,
            n_atoms=400,
            dt=0.002,
            tau=0.1,
            sigma=2.0e12,
            seed=23,
        )
        traj = _make_velocity_trajectory(velocities)

        small = GreenKubo(traj, time_step=0.002, window_size=400)
        large = GreenKubo(traj, time_step=0.002, window_size=800)
        small.run()
        large.run()

        assert small.diffusion_coefficient == pytest.approx(
            large.diffusion_coefficient,
            rel=0.05,
        )

    def test_subsampling_dt_invariance(self):
        """
        Sub-sampling the trajectory by a factor of two while doubling
        the time step reproduces the diffusion coefficient, guarding
        the ps time-step handling of the running integral.
        """
        dt = 0.001
        velocities = _ou_velocities(
            n_frames=4000,
            n_atoms=400,
            dt=dt,
            tau=0.1,
            sigma=2.0e12,
            seed=29,
        )

        full = GreenKubo(
            _make_velocity_trajectory(velocities),
            time_step=dt,
            window_size=800,
        )
        subsampled = GreenKubo(
            _make_velocity_trajectory(velocities[::2]),
            time_step=2.0 * dt,
            window_size=400,
        )
        full.run()
        subsampled.run()

        assert full.diffusion_coefficient == pytest.approx(
            subsampled.diffusion_coefficient,
            rel=0.05,
        )
