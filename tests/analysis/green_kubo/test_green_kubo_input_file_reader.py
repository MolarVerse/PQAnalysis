"""
Tests for the GreenKuboInputFileReader class.
"""

import pytest

from PQAnalysis.analysis.green_kubo import GreenKuboInputFileReader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError

from .. import pytestmark  # pylint: disable=unused-import

# pylint: disable=protected-access



def _write_input_file(filename, content):
    """
    Writes an input file with the given content.
    """
    with open(filename, "w", encoding="utf-8") as file:
        file.write(content)



MINIMAL_INPUT = (
    "traj_files = traj.vel\n"
    "target_selection = Ar\n"
    "out_file = green_kubo.dat\n"
    "time_step = 0.002\n"
)



class TestGreenKuboInputFileReader:

    """
    Tests for the GreenKuboInputFileReader class.
    """

    def test_minimal_input(self, tmpdir):  # pylint: disable=unused-argument
        """
        A minimal input file provides the required keys and None for
        all optional keys.
        """
        _write_input_file("gk.in", MINIMAL_INPUT)

        reader = GreenKuboInputFileReader("gk.in")
        reader.read()

        assert reader.traj_files == ["traj.vel"]
        assert reader.target_selection == "Ar"
        assert reader.out_file == "green_kubo.dat"
        assert reader.time_step == 0.002

        assert reader.window is None
        assert reader.gap is None
        assert reader.method is None
        assert reader.coefficient is None
        assert reader.fit_start is None
        assert reader.fit_stop is None
        assert reader.n_blocks is None
        assert reader.log_file is None

    def test_full_input(self, tmpdir):  # pylint: disable=unused-argument
        """
        All optional keys are parsed.
        """
        _write_input_file(
            "gk.in",
            (
                "traj_files = [traj_1.vel, traj_2.vel]\n"
                "target_selection = all\n"
                "out_file = gk.dat\n"
                "time_step = 0.002\n"
                "window = 500\n"
                "gap = 3\n"
                "method = direct\n"
                "coefficient = diffusion\n"
                "fit_start = 0.4\n"
                "fit_stop = 0.9\n"
                "n_blocks = 8\n"
                "log_file = gk.log\n"
                "use_full_atom_info = True\n"
            ),
        )

        reader = GreenKuboInputFileReader("gk.in")
        reader.read()

        assert reader.traj_files == ["traj_1.vel", "traj_2.vel"]
        assert reader.target_selection == "all"
        assert reader.out_file == "gk.dat"
        assert reader.time_step == 0.002
        assert reader.window == 500
        assert reader.gap == 3
        assert reader.method == "direct"
        assert reader.coefficient == "diffusion"
        assert reader.fit_start == 0.4
        assert reader.fit_stop == 0.9
        assert reader.n_blocks == 8
        assert reader.log_file == "gk.log"
        assert reader.use_full_atom_info is True

    def test_missing_required_key(self, tmpdir):  # pylint: disable=unused-argument
        """
        A missing required key (time_step) raises an InputFileError.
        """
        _write_input_file(
            "gk.in",
            (
                "traj_files = traj.vel\n"
                "target_selection = Ar\n"
                "out_file = green_kubo.dat\n"
            ),
        )

        reader = GreenKuboInputFileReader("gk.in")

        with pytest.raises(InputFileError):
            reader.read()

    def test_documentation_appended(self):
        """
        The input keys documentation is appended to the class
        docstring.
        """
        assert "Green-Kubo analysis input file" in \
            GreenKuboInputFileReader.__doc__
