"""
Tests for the MSDDataWriter and MSDLogWriter classes.
"""

import numpy as np

from PQAnalysis.analysis.msd import MSD, MSDDataWriter, MSDLogWriter
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell
from PQAnalysis.traj import Trajectory
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception



def _make_msd(time_step=None):
    """
    Builds a small MSD analysis object for writer tests.
    """
    rng = np.random.default_rng(42)
    cell = Cell(10.0, 10.0, 10.0)
    atoms = [Atom("O"), Atom("H")]

    systems = [
        AtomicSystem(
            atoms=atoms,
            pos=rng.uniform(0.0, 10.0, (2, 3)),
            cell=cell
        ) for _ in range(20)
    ]

    return MSD(
        Trajectory(systems),
        "O",
        window=5,
        gap=5,
        time_step=time_step
    )



class TestMSDDataWriter:

    """
    Tests for the MSDDataWriter class.
    """

    def test__type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("filename", 1.0, str),
            exception=PQTypeError,
            function=MSDDataWriter,
            filename=1.0,
        )

    def test_write_legacy_format(self, tmp_path):
        data = (
            np.array([0, 1]),
            np.array([0.0, 0.08081983]),
            np.array([0.0, 0.09092513]),
            np.array([0.0, 0.09018468]),
            np.array([0.0, 0.26192964]),
        )

        out_file = tmp_path / "msd.dat"
        writer = MSDDataWriter(str(out_file))

        writer.write(data)

        # exactly the legacy Diffcalc output format
        assert out_file.read_text(encoding="utf-8") == (
            "       0      0.00000000     0.00000000     0.00000000\n"
            "       1      0.08081983     0.09092513     0.09018468\n"
        )



class TestMSDLogWriter:

    """
    Tests for the MSDLogWriter class.
    """

    def test__type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "filename",
                1.0,
                str | None
            ),
            exception=PQTypeError,
            function=MSDLogWriter,
            filename=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("msd", 1.0, MSD),
            exception=PQTypeError,
            function=MSDLogWriter("test.out").write_before_run,
            msd=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("msd", 1.0, MSD),
            exception=PQTypeError,
            function=MSDLogWriter("test.out").write_after_run,
            msd=1.0,
        )

    def test_write_before_run(self, tmp_path):
        msd = _make_msd()
        log_file = tmp_path / "msd.log"
        writer = MSDLogWriter(str(log_file))

        writer.write_before_run(msd)

        contents = log_file.read_text(encoding="utf-8")

        assert "MSD calculation:" in contents
        assert "Window size (frames): 5" in contents
        assert "Origin gap (frames):  5" in contents
        assert "Number of origins:    3" in contents
        assert "Number of frames: 20" in contents
        assert "Number of atoms:  2" in contents
        assert "Target selection:" in contents
        assert "total number of atoms in target selection: 1" in contents
        assert "Time step" not in contents

    def test_write_before_run_with_time_step(self, tmp_path):
        msd = _make_msd(time_step=0.5)
        log_file = tmp_path / "msd.log"
        writer = MSDLogWriter(str(log_file))

        writer.write_before_run(msd)

        contents = log_file.read_text(encoding="utf-8")

        assert "Time step:  0.5 ps" in contents
        assert "Fit window: last 2 points" in contents

    def test_write_after_run(self, tmp_path):
        msd = _make_msd()
        msd.run()

        log_file = tmp_path / "msd.log"
        writer = MSDLogWriter(str(log_file))

        writer.write_after_run(msd)

        contents = log_file.read_text(encoding="utf-8")

        assert "Diffusion coefficients" not in contents
        assert "Elapsed time:" in contents

    def test_write_after_run_with_fit_results(self, tmp_path):
        msd = _make_msd(time_step=0.5)
        msd.run()

        log_file = tmp_path / "msd.log"
        writer = MSDLogWriter(str(log_file))

        writer.write_after_run(msd)

        contents = log_file.read_text(encoding="utf-8")

        assert "Diffusion coefficients (Einstein relation):" in contents
        assert "D_x" in contents
        assert "D_y" in contents
        assert "D_z" in contents
        assert "D_total" in contents
        assert "m^2/s" in contents
        assert "R^2" in contents
        assert "Elapsed time:" in contents
