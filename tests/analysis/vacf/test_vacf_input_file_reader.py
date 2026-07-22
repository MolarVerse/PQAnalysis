"""
Tests for the VACFInputFileReader class.
"""

import pytest

from PQAnalysis.analysis.vacf import VACFInputFileReader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception

# pylint: disable=protected-access



def _write_input_file(filename, content):
    """
    Writes an input file with the given content.
    """
    with open(filename, "w", encoding="utf-8") as file:
        file.write(content)



MINIMAL_INPUT = (
    "traj_files = traj.vel\n"
    "target_selection = O\n"
    "out_file = vacf.dat\n"
    "time_step = 0.002\n"
)



class TestVACFInputFileReader:

    """
    Tests for the VACFInputFileReader class.
    """

    def test_minimal_input(self, tmpdir):  # pylint: disable=unused-argument
        """
        A minimal input file provides the required keys and None for
        all optional keys.
        """
        _write_input_file("vacf.in", MINIMAL_INPUT)

        reader = VACFInputFileReader("vacf.in")
        reader.read()

        assert reader.traj_files == ["traj.vel"]
        assert reader.target_selection == "O"
        assert reader.out_file == "vacf.dat"
        assert reader.time_step == 0.002

        assert reader.window is None
        assert reader.gap is None
        assert reader.method is None
        assert reader.charge_file is None
        assert reader.charge_files is None
        assert reader.spectrum_file is None
        assert reader.ftsize is None
        assert reader.window_function is None
        assert reader.window_param is None
        assert reader.window_start is None
        assert reader.window_stop is None
        assert reader.windowed_out_file is None
        assert reader.log_file is None

    def test_full_input(self, tmpdir):  # pylint: disable=unused-argument
        """
        All optional keys are parsed.
        """
        _write_input_file(
            "vacf.in",
            (
                "traj_files = [traj_1.vel, traj_2.vel]\n"
                "target_selection = all\n"
                "out_file = vacf.dat\n"
                "time_step = 0.002\n"
                "window = 100\n"
                "gap = 5\n"
                "method = fft\n"
                "charge_file = charges.dat\n"
                "spectrum_file = spectrum.dat\n"
                "ftsize = 256\n"
                "window_function = blackman\n"
                "window_param = 20.0\n"
                "window_start = 0.02\n"
                "window_stop = 0.15\n"
                "windowed_out_file = windowed.dat\n"
                "log_file = vacf.log\n"
                "use_full_atom_info = True\n"
            ),
        )

        reader = VACFInputFileReader("vacf.in")
        reader.read()

        assert reader.traj_files == ["traj_1.vel", "traj_2.vel"]
        assert reader.target_selection == "all"
        assert reader.out_file == "vacf.dat"
        assert reader.time_step == 0.002
        assert reader.window == 100
        assert reader.gap == 5
        assert reader.method == "fft"
        assert reader.charge_file == "charges.dat"
        assert reader.spectrum_file == "spectrum.dat"
        assert reader.ftsize == 256
        assert reader.window_function == "blackman"
        assert reader.window_param == 20.0
        assert reader.window_start == 0.02
        assert reader.window_stop == 0.15
        assert reader.windowed_out_file == "windowed.dat"
        assert reader.log_file == "vacf.log"
        assert reader.use_full_atom_info is True

    def test_charge_files_input(self, tmpdir):  # pylint: disable=unused-argument
        """
        Charge trajectory files are parsed as a file list.
        """
        _write_input_file(
            "vacf.in",
            MINIMAL_INPUT + "charge_files = [traj_1.chrg, traj_2.chrg]\n",
        )

        reader = VACFInputFileReader("vacf.in")
        reader.read()

        assert reader.charge_files == ["traj_1.chrg", "traj_2.chrg"]

    def test_missing_required_key(self, tmpdir):  # pylint: disable=unused-argument
        """
        A missing required key raises an InputFileError.
        """
        _write_input_file(
            "vacf.in",
            (
                "traj_files = traj.vel\n"
                "target_selection = O\n"
                "out_file = vacf.dat\n"
            ),
        )

        reader = VACFInputFileReader("vacf.in")

        with pytest.raises(InputFileError):
            reader.read()

    def test_both_charge_sources(self, tmpdir, caplog):  # pylint: disable=unused-argument
        """
        A static charge file and charge trajectory files cannot be
        combined.
        """
        _write_input_file(
            "vacf.in",
            (
                MINIMAL_INPUT + "charge_file = charges.dat\n"
                "charge_files = [traj_1.chrg]\n"
            ),
        )

        reader = VACFInputFileReader("vacf.in")

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACFInputFileReader",
            logging_level="ERROR",
            message_to_test=(
                "The keys 'charge_file' and 'charge_files' cannot be "
                "used at the same time. Please provide only one charge "
                "source for the charge-flux mode."
            ),
            exception=InputFileError,
            function=reader.read,
        )

    def test_windowed_out_file_without_spectrum_file(self, tmpdir, caplog):  # pylint: disable=unused-argument
        """
        A windowed output file requires a spectrum file.
        """
        _write_input_file(
            "vacf.in",
            MINIMAL_INPUT + "windowed_out_file = windowed.dat\n",
        )

        reader = VACFInputFileReader("vacf.in")

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACFInputFileReader",
            logging_level="ERROR",
            message_to_test=(
                "The key 'windowed_out_file' can only be used together "
                "with the key 'spectrum_file'."
            ),
            exception=InputFileError,
            function=reader.read,
        )

    def test_unknown_window_function(self, tmpdir, caplog):  # pylint: disable=unused-argument
        """
        An unknown window function raises an InputFileError.
        """
        _write_input_file(
            "vacf.in",
            MINIMAL_INPUT + "window_function = hamming\n",
        )

        reader = VACFInputFileReader("vacf.in")

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="VACFInputFileReader",
            logging_level="ERROR",
            message_to_test=(
                "Unknown window function 'hamming'. Possible window "
                "functions are: none, exponential, hann, blackman."
            ),
            exception=InputFileError,
            function=reader.read,
        )
