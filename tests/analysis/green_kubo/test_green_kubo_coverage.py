"""
Additional coverage tests for the GreenKubo class that exercise the
remaining rarely-taken branches:

* the non-raw ``TrajectoryReader`` dispatch (a non-VEL reader that still
  provides velocities is read via ``frame_generator``), and
* the raw fast-path atom-count mismatch guard.
"""

import numpy as np

from PQAnalysis.analysis.green_kubo import GreenKubo
from PQAnalysis.analysis.green_kubo.exceptions import GreenKuboError
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.io.traj_file.trajectory_writer import TrajectoryWriter
from PQAnalysis.traj import Trajectory, TrajectoryFormat

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception

# pylint: disable=protected-access



def _write_vel_file(filename, n_atoms, n_frames, value=1.0):
    """
    Writes an xyz-style velocity file with a constant per-frame atom
    count.
    """
    with open(filename, "w", encoding="utf-8") as file:
        for _ in range(n_frames):
            file.write(f"{n_atoms} 10.0 10.0 10.0\n")
            file.write("comment\n")
            for _ in range(n_atoms):
                file.write(f"Ar {value} {value} {value}\n")



class TestNonRawTrajectoryReaderDispatch:

    """
    Tests for the non-raw TrajectoryReader dispatch of the velocity
    trajectory (extended xyz readers bypass the raw fast path).
    """

    def test_velocity_trajectory_reader_non_raw(self, tmp_path):
        """
        A velocity TrajectoryReader whose format is not ``VEL`` (here an
        extended xyz reader that still provides velocities) is read via
        ``calculate_number_of_frames_per_file`` and ``frame_generator``
        instead of the raw fast path.
        """
        rng = np.random.default_rng(11)
        velocities = rng.standard_normal((30, 3, 3)) * 1.0e12

        atoms = [Atom(name) for name in ("O", "H", "H")]
        systems = [
            AtomicSystem(atoms=atoms, pos=np.zeros((3, 3)), vel=frame)
            for frame in velocities
        ]
        vel_file = str(tmp_path / "vel.extxyz")
        TrajectoryWriter(vel_file).write(
            Trajectory(systems),
            traj_type=TrajectoryFormat.EXTXYZ,
        )

        reader = TrajectoryReader(vel_file)
        assert reader.traj_format != TrajectoryFormat.VEL

        analysis = GreenKubo(reader, time_step=0.1, window_size=8)

        assert analysis._raw_reader is None
        assert analysis._frame_generator is not None
        assert analysis.n_frames == 30

        lag_times, cvv, d_running = analysis.run()

        assert len(cvv) == 9
        assert cvv[0] > 0.0
        assert d_running[0] == 0.0
        assert np.allclose(lag_times, np.arange(9) * 0.1, rtol=1e-14)



class TestRawFrameAtomMismatch:

    """
    Tests for the raw fast-path atom-count mismatch guard.
    """

    def test_raw_frame_atom_count_mismatch(self, tmp_path, caplog):
        """
        A raw velocity frame that provides a different number of atoms
        than the first frame raises a GreenKuboError. The two split
        files keep the cheap frame counting consistent while the second
        file introduces the mismatch.
        """
        first_file = str(tmp_path / "first.vel")
        second_file = str(tmp_path / "second.vel")
        _write_vel_file(first_file, n_atoms=1, n_frames=15)
        _write_vel_file(second_file, n_atoms=2, n_frames=5)

        reader = TrajectoryReader(
            [first_file, second_file],
            traj_format=TrajectoryFormat.VEL,
        )

        analysis = GreenKubo(reader, time_step=0.1, window_size=10)

        assert analysis._raw_reader is not None
        assert analysis.n_atoms == 1

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
            function=analysis.run,
        )
