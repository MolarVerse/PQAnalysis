"""
Tests for the Green-Kubo API function.
"""

import os

from pathlib import Path

import numpy as np
import pytest

from PQAnalysis.analysis.green_kubo.api import green_kubo
from PQAnalysis.io.exceptions import FileWritingModeError
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception



def _write_velocity_file(path, n_frames=200, n_atoms=8, seed=17):
    """
    Writes a small Ornstein-Uhlenbeck velocity trajectory (.vel) in
    Angstrom / s.
    """
    rng = np.random.default_rng(seed)
    phi = np.exp(-0.05)
    noise_scale = 2.0 * np.sqrt(1.0 - phi * phi)

    vel = np.zeros((n_frames, n_atoms, 3))
    vel[0] = 2.0 * rng.standard_normal((n_atoms, 3))
    for n in range(1, n_frames):
        vel[n] = phi * vel[n - 1] + noise_scale * rng.standard_normal(
            (n_atoms, 3)
        )

    vel *= 1.0e12

    lines = []
    for frame in vel:
        lines.append(f"{n_atoms} 100.0 100.0 100.0\n\n")
        for x, y, z in frame:
            lines.append(f"Ar {x:.8e} {y:.8e} {z:.8e}\n")

    with open(path, "w", encoding="utf-8") as file:
        file.writelines(lines)



def _write_file(filename, content):
    """
    Writes a text file with the given content.
    """
    with open(filename, "w", encoding="utf-8") as file:
        file.write(content)



class TestGreenKuboAPI:

    """
    Tests for the input file based green_kubo api function.
    """

    def test_wrong_param_types(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("input_file", 1, str),
            exception=PQTypeError,
            function=green_kubo,
            input_file=1,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "md_format",
                1,
                "PQAnalysis.traj.formats.MDEngineFormat | str",
            ),
            exception=PQTypeError,
            function=green_kubo,
            input_file="test",
            md_format=1,
        )

    def test_green_kubo_run(self, tmpdir):  # pylint: disable=unused-argument
        """
        The api runs end-to-end from a velocity file and writes the
        running-integral data and the plateau diffusion coefficient.
        """
        _write_velocity_file("traj.vel")
        _write_file(
            "gk.in",
            (
                "traj_files = traj.vel\n"
                "target_selection = Ar\n"
                "out_file = gk.dat\n"
                "log_file = gk.log\n"
                "time_step = 0.005\n"
                "window = 60\n"
                "n_blocks = 3\n"
            ),
        )

        green_kubo("gk.in")

        out_file = Path("gk.dat")
        log_file = Path("gk.log")

        assert out_file.is_file()
        assert log_file.is_file()

        data = np.loadtxt("gk.dat")
        assert data.shape == (61, 3)
        # lag-time axis in ps
        assert np.allclose(data[:, 0], np.arange(61) * 0.005)
        # running integral starts at zero
        assert data[0, 2] == 0.0

        log_contents = log_file.read_text(encoding="utf-8")
        assert "Green-Kubo self-diffusion coefficient calculation:" \
            in log_contents
        assert "m^2/s" in log_contents
        assert "cm^2/s" in log_contents
        assert "(n_blocks=3)" in log_contents
        assert "Elapsed time:" in log_contents

    def test_existing_out_file_fails_before_run(self, tmpdir):  # pylint: disable=unused-argument
        """
        A pre-existing output file raises a FileWritingModeError before
        the analysis runs, so that no computation is lost.
        """
        _write_velocity_file("traj.vel")
        _write_file(
            "gk.in",
            (
                "traj_files = traj.vel\n"
                "target_selection = Ar\n"
                "out_file = gk.dat\n"
                "time_step = 0.005\n"
                "window = 60\n"
            ),
        )
        _write_file("gk.dat", "already there\n")

        with pytest.raises(FileWritingModeError):
            green_kubo("gk.in")

        assert not os.path.exists("gk.log")
