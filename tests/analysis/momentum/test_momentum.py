"""
Tests for the Momentum analysis class.
"""

import sys

import numpy as np
import pytest

from PQAnalysis import config
from PQAnalysis.analysis.momentum import Momentum
from PQAnalysis.analysis.momentum.exceptions import MomentumError
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.traj import Trajectory

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception

# pylint: disable=protected-access

momentum_module = sys.modules[Momentum.__module__]

# hand-computed |P| references for tests/data/momentum/two_frames.vel
# (velocities are exactly representable in float32, masses
# O = 15.9994 amu and H = 1.00794 amu):
# frame 1: P = 15.9994*(0.5, -0.25, 1.0) + 1.00794*(2.0, 1.5, -0.5)
#              + 1.00794*(-1.0, 0.25, 0.75)
# frame 2: P = 15.9994*(1.0, 0.0, 0.0) + 1.00794*(0.0, 2.0, 0.0)
#              + 1.00794*(0.0, 0.0, 4.0)
TWO_FRAMES_NORMS = [18.714822669473786, 16.622264022448928]
TWO_FRAMES_NORMS_OXYGEN = [18.329615393469116, 15.9994]

# masses used by the independent gas3.vel reference parser below
# (identical to the PQAnalysis element data for c and o)
GAS3_MASSES = {"c": 12.0107, "o": 15.9994}


def _independent_gas3_norms(filename, scale=1e-15):
    """
    Compute the scaled total momentum norms of a QMCFC velocity
    trajectory independently of the TrajectoryReader/Momentum
    pipeline.

    This minimal parser only shares the documented number semantics
    with the pipeline: every velocity component is parsed from its
    token in single precision (float32 quantization) and the
    mass-weighted sum is accumulated in float64. Any behavioral
    change in TrajectoryReader parsing or in the Momentum
    accumulation therefore makes the pipeline diverge from this
    reference.
    """
    with open(filename, encoding="utf-8") as file:
        lines = file.readlines()

    norms = []
    line_index = 0

    while line_index < len(lines):
        n_atoms = int(lines[line_index].split()[0])
        atom_lines = lines[line_index + 2:line_index + 2 + n_atoms]
        line_index += 2 + n_atoms

        momentum = np.zeros(3, dtype=np.float64)

        for atom_line in atom_lines:
            tokens = atom_line.split()
            name = tokens[0].lower()

            if name == "x":  # QMCFC dummy atom is stripped
                continue

            velocity = np.array(
                tokens[1:4], dtype=np.float32
            ).astype(np.float64)
            momentum += GAS3_MASSES[name] * velocity

        norms.append(float(np.linalg.norm(momentum)) * scale)

    return np.array(norms)



class TestMomentum:

    """
    Tests for the Momentum class.
    """

    @pytest.mark.parametrize("example_dir", ["momentum"], indirect=False)
    def test_two_frames_reference(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        The momentum norms of a hand-computed two-frame QMCFC velocity
        trajectory are reproduced. The leading dummy 'X' atom is
        stripped by the QMCFC engine format.
        """
        reader = TrajectoryReader("two_frames.vel", md_format="qmcfc")
        momentum_norms = Momentum(reader).run()

        assert np.allclose(
            momentum_norms,
            np.array(TWO_FRAMES_NORMS) * 1e-15,
            rtol=1e-12,
        )

    @pytest.mark.parametrize("example_dir", ["momentum"], indirect=False)
    def test_gas_noise_floor_reference(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        The momentum norms of the first three frames of the legacy
        equipartition.jl gas phase test trajectory (gas3.vel, the
        momentum-conserving 150-atom QMCFC trajectory gas_md_01.vel)
        match an independent single-precision reference computation.

        The total momentum of this data is zero up to floating point
        noise (single m_i * v_i terms are of the order 5e13
        amu*Angstrom/s), so the norm is pure summation noise of the
        parsing precision: PQAnalysis parses the velocities as
        float32 and reports norms of the order 5e-8 (scaled), while
        the legacy Julia tool parses as float64 and prints
        2.6407665688752317e-16, 3.625033674412557e-16 and
        2.827919968396604e-16 for these frames. The expected values
        are recomputed by _independent_gas3_norms, a minimal parser
        written independently of the TrajectoryReader/Momentum
        pipeline, so this test fails if either changes behavior.
        """
        expected_norms = _independent_gas3_norms("gas3.vel")

        reader = TrajectoryReader("gas3.vel", md_format="qmcfc")
        momentum_norms = Momentum(reader).run()

        # momentum conservation: the norms sit at the float32
        # parsing noise floor, orders of magnitude below a single
        # scaled m_i * |v_i| term (~5e-2 after the 1e-15 scaling)
        assert np.all(momentum_norms < 1e-6)

        assert np.allclose(momentum_norms, expected_norms, rtol=1e-6)

    @pytest.mark.parametrize("example_dir", ["momentum"], indirect=False)
    def test_selection(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        With a selection only the selected atoms contribute to the
        total momentum.
        """
        reader = TrajectoryReader("two_frames.vel", md_format="qmcfc")
        momentum_norms = Momentum(reader, selection="O").run()

        assert np.allclose(
            momentum_norms,
            np.array(TWO_FRAMES_NORMS_OXYGEN) * 1e-15,
            rtol=1e-12,
        )

    @pytest.mark.parametrize("example_dir", ["momentum"], indirect=False)
    def test_scale(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        A custom scaling factor replaces the default 1e-15 conversion.
        """
        reader = TrajectoryReader("two_frames.vel", md_format="qmcfc")
        momentum_norms = Momentum(reader, scale=1.0).run()

        assert np.allclose(
            momentum_norms, TWO_FRAMES_NORMS, rtol=1e-12
        )

    def test_zero_drift(self):
        """
        Symmetric velocities of identical atoms cancel to a total
        momentum of exactly zero.
        """
        system = AtomicSystem(
            atoms=[Atom("H"), Atom("H")],
            vel=np.array([[1.5, -2.0, 3.25], [-1.5, 2.0, -3.25]]),
        )
        momentum_norms = Momentum(Trajectory([system])).run()

        assert momentum_norms.shape == (1, )
        assert momentum_norms[0] == 0.0

    def test_n_frames(self):
        """
        The number of analyzed frames is exposed after run().
        """
        system = AtomicSystem(
            atoms=[Atom("H")], vel=np.array([[1.0, 0.0, 0.0]])
        )
        momentum = Momentum(Trajectory([system, system]))

        assert momentum.n_frames == 0

        momentum.run()

        assert momentum.n_frames == 2

    def test_empty_trajectory(self, caplog):
        """
        An empty trajectory is rejected.
        """
        assert_logging_with_exception(
            caplog=caplog,
            logging_name=Momentum.__qualname__,
            logging_level="ERROR",
            message_to_test="Trajectory cannot be of length 0.",
            exception=MomentumError,
            function=Momentum,
            traj=Trajectory(),
        )

    def test_empty_selection(self, caplog):
        """
        A selection that does not select any atoms is rejected
        instead of silently producing all-zero momentum norms.
        """
        system = AtomicSystem(
            atoms=[Atom("H"), Atom("H")],
            vel=np.array([[1.5, -2.0, 3.25], [-1.5, 2.0, -3.25]]),
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=Momentum.__qualname__,
            logging_level="ERROR",
            message_to_test="The selection does not select any atoms.",
            exception=MomentumError,
            function=Momentum,
            traj=Trajectory([system]),
            selection="O",
        )

    def test_unknown_mass(self, caplog):
        """
        A selection containing an atom with unknown mass is rejected.
        """
        system = AtomicSystem(
            atoms=[Atom("Foo123", use_guess_element=False)],
            vel=np.array([[1.0, 0.0, 0.0]]),
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=Momentum.__qualname__,
            logging_level="ERROR",
            message_to_test=(
                "The mass of at least one selected atom is unknown. "
                "The total momentum cannot be calculated."
            ),
            exception=MomentumError,
            function=Momentum,
            traj=Trajectory([system]),
        )

    @pytest.mark.parametrize("example_dir", ["momentum"], indirect=False)
    def test_raw_fast_path_matches_in_memory_path(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        A velocity TrajectoryReader dispatches to the raw fast-path
        stream, which produces bit-identical momentum norms compared
        to the in-memory AtomicSystem based stream.
        """
        momentum = Momentum(
            TrajectoryReader("gas3.vel", md_format="qmcfc")
        )

        assert momentum._raw_reader is not None

        raw_norms = momentum.run()

        traj = TrajectoryReader("gas3.vel", md_format="qmcfc").read()
        in_memory = Momentum(traj)

        assert in_memory._raw_reader is None

        in_memory_norms = in_memory.run()

        assert np.array_equal(raw_norms, in_memory_norms)

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

        monkeypatch.setattr(momentum_module, "tqdm", fake_tqdm)

        system = AtomicSystem(
            atoms=[Atom("H")], vel=np.array([[1.0, 0.0, 0.0]])
        )

        monkeypatch.setattr(config, "with_progress_bar", False)
        Momentum(Trajectory([system, system])).run()
        assert captured["disable"] is True

        captured.clear()

        monkeypatch.setattr(config, "with_progress_bar", True)
        Momentum(Trajectory([system, system])).run()
        assert captured["disable"] is False

    def test_missing_velocities(self, caplog):
        """
        A trajectory without velocity information is rejected.
        """
        system = AtomicSystem(
            atoms=[Atom("H")], pos=np.array([[0.0, 0.0, 0.0]])
        )
        momentum = Momentum(Trajectory([system]))

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=Momentum.__qualname__,
            logging_level="ERROR",
            message_to_test=(
                "The trajectory does not contain velocity "
                "information for all atoms. Please provide a "
                "velocity trajectory."
            ),
            exception=MomentumError,
            function=momentum.run,
        )
