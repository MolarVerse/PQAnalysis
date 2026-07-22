"""
Tests for the Momentum trajectory-source dispatch branches.

These tests exercise the construction-time dispatch of the Momentum
class that is not covered by ``test_momentum.py``: the explicit
``use_full_atom_info=None`` default fallback and the lazy
TrajectoryReader branch for a non-VEL (here extended xyz) velocity
trajectory, which streams frames through ``frame_generator`` instead
of the raw VEL fast path.
"""

import numpy as np

from PQAnalysis.analysis.momentum import Momentum
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.traj import Trajectory, TrajectoryFormat

from .. import pytestmark  # pylint: disable=unused-import

# pylint: disable=protected-access

# extended xyz velocity trajectory (parsed in float64 by the reader);
# hand-computed |P| = m_H * |sum_i v_i| references with m_H = 1.00794
_EXTXYZ_VEL_TRAJECTORY = (
    "2\n"
    "Properties=species:S:1:pos:R:3:vel:R:3\n"
    "H 0.0 0.0 0.0 1.0 2.0 3.0\n"
    "H 1.0 1.0 1.0 0.5 0.5 0.5\n"
    "2\n"
    "Properties=species:S:1:pos:R:3:vel:R:3\n"
    "H 0.0 0.0 0.0 -1.0 0.0 2.0\n"
    "H 1.0 1.0 1.0 0.0 3.0 -1.0\n"
)

# per-frame summed velocity vectors of the trajectory above
_EXTXYZ_VEL_SUMS = np.array([[1.5, 2.5, 3.5], [-1.0, 3.0, 1.0]])



class TestMomentumDispatch:

    """
    Tests for the Momentum construction-time source dispatch.
    """

    def test_use_full_atom_info_none_uses_default(self):
        """
        Passing ``use_full_atom_info=None`` explicitly falls back to the
        class default (``False``) instead of storing ``None``.
        """
        system = AtomicSystem(
            atoms=[Atom("H")], vel=np.array([[1.0, 0.0, 0.0]])
        )

        momentum = Momentum(
            Trajectory([system]), use_full_atom_info=None
        )

        assert momentum.use_full_atom_info is False
        assert (
            momentum.use_full_atom_info == Momentum._use_full_atom_default
        )

    def test_non_vel_reader_uses_frame_generator(self, tmp_path):
        """
        A TrajectoryReader whose format is not VEL (an extended xyz
        velocity trajectory) does not take the raw VEL fast path but
        the lazy frame_generator branch, and still computes the correct
        scaled total momentum norms.
        """
        traj_file = tmp_path / "vel.extxyz"
        traj_file.write_text(_EXTXYZ_VEL_TRAJECTORY, encoding="utf-8")

        reader = TrajectoryReader(str(traj_file), traj_format="EXTXYZ")

        assert reader.traj_format != TrajectoryFormat.VEL

        momentum = Momentum(reader)

        # the non-VEL reader lands in the frame_generator branch, not
        # the raw VEL fast path
        assert momentum._raw_reader is None
        assert momentum.frame_generator is not None
        assert momentum._n_frames_total == 2

        norms = momentum.run()

        mass = Atom("H").mass
        expected = mass * np.linalg.norm(_EXTXYZ_VEL_SUMS, axis=1) * 1e-15

        assert np.allclose(norms, expected, rtol=1e-12)

    def test_non_vel_reader_matches_in_memory_path(self, tmp_path):
        """
        The lazy non-VEL reader branch yields bit-identical momentum
        norms to the in-memory Trajectory branch for the same frames.
        """
        traj_file = tmp_path / "vel.extxyz"
        traj_file.write_text(_EXTXYZ_VEL_TRAJECTORY, encoding="utf-8")

        reader = TrajectoryReader(str(traj_file), traj_format="EXTXYZ")
        lazy = Momentum(reader)

        assert lazy._raw_reader is None
        assert lazy.frame_generator is not None

        lazy_norms = lazy.run()

        traj = TrajectoryReader(
            str(traj_file), traj_format="EXTXYZ"
        ).read()
        in_memory = Momentum(traj)

        assert in_memory._raw_reader is None

        assert np.array_equal(lazy_norms, in_memory.run())
