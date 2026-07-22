"""
Additional coverage tests for the VACF class and the input-file based
api that exercise the remaining rarely-taken branches:

* the ``ModuleNotFoundError`` fallback import of the accumulation
  kernels,
* the non-positive ``time_step`` guard,
* the non-raw ``TrajectoryReader`` dispatch for both the velocity and
  the charge trajectory,
* the raw fast-path atom-count mismatch guard,
* the lockstep charge-stream exhaustion and shape guards,
* the fft zero-norm guard, and
* the ``window_param`` spectrum keyword forwarding of the api.
"""

import importlib
import os
import shutil
import sys

from pathlib import Path

import numpy as np

from PQAnalysis.analysis.vacf import VACF, vacf
from PQAnalysis.analysis.vacf.exceptions import VACFError
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.io.traj_file.trajectory_writer import TrajectoryWriter
from PQAnalysis.traj import Trajectory, TrajectoryFormat

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception

# pylint: disable=protected-access

#: The compiled Cython kernel of the accumulation functions.
_KERNEL = "PQAnalysis.analysis.vacf._vacf_kernel"

#: The module object of the VACF class, used for the in-process reload
#: that exercises the ModuleNotFoundError fallback import.
_vacf_module = sys.modules[VACF.__module__]



class _KernelBlocker:

    """
    A meta_path finder that hides the compiled vacf kernel so that the
    ``except ModuleNotFoundError`` fallback import is taken on reload.
    """

    def find_spec(self, name, path=None, target=None):  # pylint: disable=unused-argument
        """Raise ModuleNotFoundError for the kernel, delegate otherwise."""
        if name == _KERNEL:
            raise ModuleNotFoundError(name)



def _make_velocity_systems(velocities, names=None):
    """
    Builds a list of velocity AtomicSystems (with dummy positions).
    """
    velocities = np.asarray(velocities, dtype=float)
    n_atoms = velocities.shape[1]

    if names is None:
        names = ["O"] * n_atoms

    atoms = [Atom(name) for name in names]

    return [
        AtomicSystem(
            atoms=atoms,
            pos=np.zeros((n_atoms, 3)),
            vel=frame_velocities,
        ) for frame_velocities in velocities
    ]



def _write_vel_file(filename, n_atoms, n_frames, value=1.0):
    """
    Writes a legacy xyz-style velocity file with a constant per-frame
    atom count.
    """
    with open(filename, "w", encoding="utf-8") as file:
        for _ in range(n_frames):
            file.write(f"{n_atoms} 10.0 10.0 10.0\n")
            file.write("comment\n")
            for _ in range(n_atoms):
                file.write(f"O {value} {value} {value}\n")


#: The legacy reference data directory of the VACF analysis, resolved
#: relative to this test file so that it does not depend on the cwd.
_DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "vacf"



class TestKernelFallbackImport:

    """
    Tests for the ModuleNotFoundError fallback import of the kernels.
    """

    def test_module_not_found_fallback_import(self):
        """
        When the compiled kernel is absent, ``accumulate_frame`` and
        ``weight_frame`` are imported from the pure-python
        ``_vacf_kernel_py`` fallback.

        The vacf module is reloaded in-process with the compiled kernel
        hidden by a meta_path blocker and restored to the compiled
        version afterwards, so that the test does not affect any other
        test regardless of the execution order.
        """
        assert _vacf_module.accumulate_frame.__module__.endswith(
            "_vacf_kernel"
        )
        assert _vacf_module.weight_frame.__module__.endswith("_vacf_kernel")

        sys.modules.pop(_KERNEL, None)
        blocker = _KernelBlocker()
        sys.meta_path.insert(0, blocker)

        try:
            importlib.reload(_vacf_module)
            assert _vacf_module.accumulate_frame.__module__.endswith(
                "_vacf_kernel_py"
            )
            assert _vacf_module.weight_frame.__module__.endswith(
                "_vacf_kernel_py"
            )
        finally:
            sys.meta_path.remove(blocker)
            sys.modules.pop(_KERNEL, None)
            importlib.reload(_vacf_module)

        assert _vacf_module.accumulate_frame.__module__.endswith(
            "_vacf_kernel"
        )
        assert _vacf_module.weight_frame.__module__.endswith("_vacf_kernel")



class TestNonPositiveTimeStep:

    """
    Tests for the non-positive time step guard.
    """

    def test_zero_time_step(self, caplog):
        """
        A time step of exactly zero raises a VACFError. ``0.0`` is a
        valid ``PositiveReal`` (``>= 0.0``), so the guard is reached in
        both the debug and the release type-checking modes.
        """
        traj = Trajectory(_make_velocity_systems(np.ones((20, 1, 3))))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "The time_step must be a positive real number."
            ),
            exception=VACFError,
            function=VACF,
            traj=traj,
            window_size=5,
            time_step=0.0,
        )



class TestNonRawTrajectoryReaderDispatch:

    """
    Tests for the non-raw TrajectoryReader dispatch of the velocity and
    the charge trajectory (extended xyz readers bypass the raw
    fast paths).
    """

    def test_velocity_trajectory_reader_non_raw(self, tmp_path):
        """
        A velocity TrajectoryReader whose format is not ``VEL`` (here an
        extended xyz reader that still provides velocities) is read via
        ``calculate_number_of_frames_per_file`` and ``frame_generator``
        instead of the raw fast path.
        """
        rng = np.random.default_rng(11)
        velocities = rng.standard_normal((30, 3, 3))

        systems = _make_velocity_systems(velocities, names=["O", "H", "H"])
        vel_file = str(tmp_path / "vel.extxyz")
        TrajectoryWriter(vel_file).write(
            Trajectory(systems),
            traj_type=TrajectoryFormat.EXTXYZ,
        )

        reader = TrajectoryReader(vel_file)
        assert reader.traj_format != TrajectoryFormat.VEL

        analysis = VACF(reader, window_size=8, time_step=0.1, gap=2)

        assert analysis._raw_reader is None
        assert analysis._frame_generator is not None
        assert analysis.n_frames == 30

        time, correlation = analysis.run()

        assert correlation[0] == 1.0
        assert len(correlation) == 9
        assert np.allclose(time, np.arange(9) * 0.1, rtol=1e-14)

    def test_charge_trajectory_reader_non_raw(self, tmp_path):
        """
        A charge TrajectoryReader whose format is not ``CHARGE`` (here
        an extended xyz reader that provides charges) is read via
        ``calculate_number_of_frames_per_file`` and ``frame_generator``
        instead of the raw charge fast path.
        """
        rng = np.random.default_rng(12)
        n_frames, n_atoms = 30, 3
        velocities = rng.standard_normal((n_frames, n_atoms, 3))
        charges = 0.5 + 0.1 * rng.standard_normal((n_frames, n_atoms))

        atoms = [Atom(name) for name in ("O", "H", "H")]
        charge_systems = [
            AtomicSystem(
                atoms=atoms,
                pos=np.zeros((n_atoms, 3)),
                vel=velocities[i],
                charges=charges[i],
            ) for i in range(n_frames)
        ]
        charge_file = str(tmp_path / "chg.extxyz")
        TrajectoryWriter(charge_file).write(
            Trajectory(charge_systems),
            traj_type=TrajectoryFormat.EXTXYZ,
        )

        charge_reader = TrajectoryReader(charge_file)
        assert charge_reader.traj_format != TrajectoryFormat.CHARGE

        velocity_traj = Trajectory(
            _make_velocity_systems(velocities, names=["O", "H", "H"])
        )

        analysis = VACF(
            velocity_traj,
            window_size=8,
            time_step=0.1,
            gap=2,
            charge_traj=charge_reader,
        )

        assert analysis._raw_charge_reader is None
        assert analysis._charge_frame_generator is not None
        assert analysis.flux is True

        _, correlation = analysis.run()

        assert correlation[0] == 1.0



class TestRawFrameAtomMismatch:

    """
    Tests for the raw fast-path atom-count mismatch guard.
    """

    def test_raw_frame_atom_count_mismatch(self, tmp_path, caplog):
        """
        A raw velocity frame that provides a different number of atoms
        than the first frame raises a VACFError. The two split files
        keep the cheap frame counting consistent while the second file
        introduces the mismatch.
        """
        first_file = str(tmp_path / "first.vel")
        second_file = str(tmp_path / "second.vel")
        _write_vel_file(first_file, n_atoms=1, n_frames=15)
        _write_vel_file(second_file, n_atoms=2, n_frames=5)

        reader = TrajectoryReader(
            [first_file, second_file],
            traj_format=TrajectoryFormat.VEL,
        )

        analysis = VACF(reader, window_size=10, time_step=0.1)

        assert analysis._raw_reader is not None
        assert analysis.n_atoms == 1

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
            function=analysis.run,
        )



class TestChargeStreamGuards:

    """
    Tests for the lockstep charge-stream exhaustion and shape guards.
    """

    def _make_charge_flux_analysis(self, seed):
        """
        Builds a valid charge-flux VACF whose charge stream can be
        overridden to exercise the guards of ``_next_charges``.
        """
        rng = np.random.default_rng(seed)
        n_frames, n_atoms = 20, 3
        velocities = rng.standard_normal((n_frames, n_atoms, 3))
        charges = 0.5 + 0.1 * rng.standard_normal((n_frames, n_atoms))

        atoms = [Atom("O")] * n_atoms
        velocity_traj = Trajectory(
            _make_velocity_systems(velocities)
        )
        charge_traj = Trajectory(
            [
                AtomicSystem(atoms=atoms, charges=charges[i])
                for i in range(n_frames)
            ]
        )

        return VACF(
            velocity_traj,
            window_size=5,
            time_step=0.1,
            charge_traj=charge_traj,
        )

    def test_charge_stream_exhausted(self, caplog):
        """
        A charge stream that runs out before the velocity trajectory
        raises a VACFError.
        """
        analysis = self._make_charge_flux_analysis(seed=21)
        analysis._charge_value_stream = iter([])

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "The charge trajectory provides fewer frames than the "
                "velocity trajectory."
            ),
            exception=VACFError,
            function=analysis._next_charges,
        )

    def test_charge_frame_wrong_shape(self, caplog):
        """
        A charge frame that does not provide one charge per atom raises
        a VACFError.
        """
        analysis = self._make_charge_flux_analysis(seed=22)
        analysis._charge_value_stream = iter(
            (np.zeros(analysis.n_atoms + 1),)
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "A frame of the charge trajectory does not provide "
                "charges for all 3 atoms. Please provide a charge "
                "trajectory (e.g. .chrg files)."
            ),
            exception=VACFError,
            function=analysis._next_charges,
        )



class TestFFTZeroNorm:

    """
    Tests for the fft zero-norm guard.
    """

    def test_fft_zero_norm(self, caplog):
        """
        An all-zero velocity trajectory has a vanishing aggregate
        squared velocity norm and raises a VACFError for the fft
        estimator.
        """
        traj = Trajectory(_make_velocity_systems(np.zeros((20, 1, 3))))
        analysis = VACF(traj, window_size=5, time_step=0.1, method="fft")

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACF",
            logging_level="ERROR",
            message_to_test=(
                "The aggregate squared velocity norm of the trajectory "
                "is zero. The normalized VACF is not defined."
            ),
            exception=VACFError,
            function=analysis.run,
        )



class TestApiWindowParam:

    """
    Tests for the ``window_param`` spectrum keyword forwarding of the
    input-file based api.
    """

    def test_window_param_forwarded_to_spectrum(self, tmp_path, monkeypatch):
        """
        A ``window_param`` set in the input file is forwarded to the
        spectrum calculation and a spectrum file is produced.
        """
        shutil.copytree(_DATA_DIR, tmp_path, dirs_exist_ok=True)
        monkeypatch.chdir(tmp_path)

        with open("vacf.in", "w", encoding="utf-8") as file:
            file.write(
                "traj_files = [traj_1.vel, traj_2.vel]\n"
                "target_selection = all\n"
                "out_file = vacf_out.dat\n"
                "time_step = 0.002\n"
                "window = 100\n"
                "gap = 5\n"
                "spectrum_file = spectrum_out.dat\n"
                "window_function = exponential\n"
                "window_param = 20.0\n"
                "window_start = 0.02\n"
                "window_stop = 0.15\n"
            )

        vacf("vacf.in", md_format="qmcfc")

        assert os.path.exists("vacf_out.dat")
        assert os.path.exists("spectrum_out.dat")
