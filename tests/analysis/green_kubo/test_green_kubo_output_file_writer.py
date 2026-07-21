"""
Tests for the GreenKuboDataWriter and GreenKuboLogWriter classes.
"""

import numpy as np

from PQAnalysis.analysis.green_kubo import (
    GreenKubo,
    GreenKuboDataWriter,
    GreenKuboLogWriter,
)
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.traj import Trajectory
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception



def _make_green_kubo():
    """
    Builds a small GreenKubo analysis object for writer tests.
    """
    rng = np.random.default_rng(42)
    atoms = [Atom("Ar"), Atom("Ar")]

    systems = [
        AtomicSystem(atoms=atoms, vel=rng.standard_normal((2, 3)))
        for _ in range(40)
    ]

    return GreenKubo(
        Trajectory(systems),
        time_step=0.002,
        window_size=10,
        gap=2,
    )



class TestGreenKuboDataWriter:

    """
    Tests for the GreenKuboDataWriter class.
    """

    def test__type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("filename", 1.0, str),
            exception=PQTypeError,
            function=GreenKuboDataWriter,
            filename=1.0,
        )

    def test_write_three_columns(self, tmp_path):
        data = (
            np.array([0.0, 0.002]),
            np.array([4.0e24, 1.5e24]),
            np.array([0.0, 3.3e-9]),
        )

        out_file = tmp_path / "green_kubo.dat"
        writer = GreenKuboDataWriter(str(out_file))

        writer.write(data)

        lines = out_file.read_text(encoding="utf-8").splitlines()

        assert len(lines) == 2
        # each row has three columns: lag time, Cvv, D_running
        assert len(lines[0].split()) == 3
        assert lines[0].split()[0] == "0.000000"



class TestGreenKuboLogWriter:

    """
    Tests for the GreenKuboLogWriter class.
    """

    def test__type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "filename",
                1.0,
                str | None,
            ),
            exception=PQTypeError,
            function=GreenKuboLogWriter,
            filename=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "green_kubo",
                1.0,
                GreenKubo,
            ),
            exception=PQTypeError,
            function=GreenKuboLogWriter("test.out").write_before_run,
            green_kubo=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "green_kubo",
                1.0,
                GreenKubo,
            ),
            exception=PQTypeError,
            function=GreenKuboLogWriter("test.out").write_after_run,
            green_kubo=1.0,
        )

    def test_write_before_run(self, tmp_path):
        green_kubo = _make_green_kubo()
        log_file = tmp_path / "green_kubo.log"
        writer = GreenKuboLogWriter(str(log_file))

        writer.write_before_run(green_kubo)

        contents = log_file.read_text(encoding="utf-8")

        assert "Green-Kubo self-diffusion coefficient calculation:" in contents
        assert "Coefficient:          diffusion" in contents
        assert "Method:               fft" in contents
        assert "Window size (frames): 10" in contents
        assert "Origin gap (frames):  2" in contents
        assert "Time step:            0.002 ps" in contents
        assert "Number of frames: 40" in contents
        assert "Number of atoms:  2" in contents
        assert "Number of blocks:     5" in contents
        assert "total number of atoms in target selection: 2" in contents

    def test_write_after_run(self, tmp_path):
        green_kubo = _make_green_kubo()
        green_kubo.run()

        log_file = tmp_path / "green_kubo.log"
        writer = GreenKuboLogWriter(str(log_file))

        writer.write_after_run(green_kubo)

        contents = log_file.read_text(encoding="utf-8")

        assert "Number of origins:" in contents
        assert "Green-Kubo self-diffusion coefficient (plateau):" in contents
        assert "m^2/s" in contents
        assert "cm^2/s" in contents
        # the +/- uncertainty is the block-averaged standard error
        n_blocks_used = len(green_kubo.block_diffusion_coefficients)
        assert f"(n_blocks={n_blocks_used})" in contents
        assert "Plateau spread" in contents
        assert "Elapsed time:" in contents

    def test_log_to_stdout(self, capsys):
        green_kubo = _make_green_kubo()
        green_kubo.run()

        writer = GreenKuboLogWriter(None)
        writer.write_before_run(green_kubo)
        writer.write_after_run(green_kubo)

        captured = capsys.readouterr()

        assert "Green-Kubo self-diffusion coefficient" in captured.out
