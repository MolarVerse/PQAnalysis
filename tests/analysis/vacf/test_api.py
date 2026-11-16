"""
Tests for the API functions of the VACF analysis.
"""

import importlib
import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from PQAnalysis.analysis.vacf import read_static_charges, vacf
from PQAnalysis.analysis.vacf.exceptions import VACFError
from PQAnalysis.io.exceptions import FileWritingModeError

from .. import pytestmark  # pylint: disable=unused-import

EXAMPLES_WATER = Path(__file__).resolve().parents[3] / "examples" / "water"

# pylint: disable=protected-access



def _write_file(filename, content):
    """
    Writes a file with the given content.
    """
    with open(filename, "w", encoding="utf-8") as file:
        file.write(content)



class TestReadStaticCharges:

    """
    Tests for the read_static_charges function.
    """

    def test_read_pq_format(self, tmpdir):  # pylint: disable=unused-argument
        """
        For the PQ format all charge entries are read.
        """
        _write_file(
            "charges.dat",
            (
                "3\n"
                "# comment\n"
                "O -0.82\n"
                "H 0.41\n"
                "H 0.41\n"
            ),
        )

        charges = read_static_charges("charges.dat")

        assert np.allclose(charges, [-0.82, 0.41, 0.41], atol=1e-14)

    def test_read_qmcfc_format_strips_dummy_atom(self, tmpdir):  # pylint: disable=unused-argument
        """
        For the QMCFC format the leading dummy 'X' entry is stripped.
        """
        _write_file(
            "charges.dat",
            (
                "4\n"
                "# comment\n"
                "X 0.0\n"
                "O -0.82\n"
                "H 0.41\n"
                "H 0.41\n"
            ),
        )

        charges = read_static_charges("charges.dat", md_format="qmcfc")

        assert np.allclose(charges, [-0.82, 0.41, 0.41], atol=1e-14)

    def test_qmcfc_format_without_dummy_atom(self, tmpdir):  # pylint: disable=unused-argument
        """
        A QMCFC charge file without a leading dummy 'X' entry raises a
        VACFError.
        """
        _write_file(
            "charges.dat",
            (
                "2\n"
                "# comment\n"
                "O -0.82\n"
                "H 0.41\n"
            ),
        )

        with pytest.raises(VACFError, match="dummy 'X' atom"):
            read_static_charges("charges.dat", md_format="qmcfc")

    def test_qmcfc_format_zero_atoms(self, tmpdir):  # pylint: disable=unused-argument
        """
        A QMCFC charge file declaring zero atoms raises a VACFError
        instead of a raw IndexError (it cannot contain the dummy 'X'
        atom).
        """
        _write_file(
            "charges.dat",
            (
                "0\n"
                "# comment\n"
            ),
        )

        with pytest.raises(VACFError, match="dummy 'X' atom"):
            read_static_charges("charges.dat", md_format="qmcfc")

    def test_invalid_header(self, tmpdir):  # pylint: disable=unused-argument
        """
        An unparsable header raises a VACFError.
        """
        _write_file("charges.dat", "not_a_number\n# comment\n")

        with pytest.raises(VACFError, match="number of atoms"):
            read_static_charges("charges.dat")

    def test_too_few_charge_lines(self, tmpdir):  # pylint: disable=unused-argument
        """
        Fewer charge lines than announced in the header raise a
        VACFError.
        """
        _write_file(
            "charges.dat",
            (
                "3\n"
                "# comment\n"
                "O -0.82\n"
                "H 0.41\n"
            ),
        )

        with pytest.raises(VACFError, match="provides only"):
            read_static_charges("charges.dat")

    def test_invalid_charge_line(self, tmpdir):  # pylint: disable=unused-argument
        """
        An unparsable charge line raises a VACFError.
        """
        _write_file(
            "charges.dat",
            (
                "2\n"
                "# comment\n"
                "O -0.82\n"
                "H not_a_number\n"
            ),
        )

        with pytest.raises(VACFError, match="Could not parse"):
            read_static_charges("charges.dat")



class TestVACFApi:

    """
    Tests for the input file based vacf api function.
    """

    @pytest.mark.parametrize("example_dir", ["vacf"], indirect=False)
    def test_existing_spectrum_file_fails_before_run(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        A pre-existing spectrum output file raises a
        FileWritingModeError before the VACF analysis is run, so that
        no computation is lost: the VACF output file must not have
        been written yet.
        """
        _write_file(
            "vacf.in",
            (
                "traj_files = [traj_1.vel, traj_2.vel]\n"
                "target_selection = all\n"
                "out_file = vacf_out.dat\n"
                "time_step = 0.002\n"
                "window = 100\n"
                "gap = 5\n"
                "spectrum_file = spectrum_out.dat\n"
                "windowed_out_file = windowed_out.dat\n"
                "window_function = hann\n"
            ),
        )

        _write_file("spectrum_out.dat", "already there\n")

        with pytest.raises(FileWritingModeError):
            vacf("vacf.in", md_format="qmcfc")

        assert not os.path.exists("vacf_out.dat")
        assert not os.path.exists("windowed_out.dat")

    @pytest.mark.parametrize("example_dir", ["vacf"], indirect=False)
    def test_existing_windowed_out_file_fails_before_run(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        A pre-existing windowed output file raises a
        FileWritingModeError before the VACF analysis is run.
        """
        _write_file(
            "vacf.in",
            (
                "traj_files = [traj_1.vel, traj_2.vel]\n"
                "target_selection = all\n"
                "out_file = vacf_out.dat\n"
                "time_step = 0.002\n"
                "window = 100\n"
                "gap = 5\n"
                "spectrum_file = spectrum_out.dat\n"
                "windowed_out_file = windowed_out.dat\n"
                "window_function = hann\n"
            ),
        )

        _write_file("windowed_out.dat", "already there\n")

        with pytest.raises(FileWritingModeError):
            vacf("vacf.in", md_format="qmcfc")

        assert not os.path.exists("vacf_out.dat")
        assert not os.path.exists("spectrum_out.dat")



def _fail(*_args, **_kwargs):
    raise AssertionError("trajectory was read")



class TestVACFAPIExistingOutput:

    """
    Tests that existing output files raise before the trajectory is read.
    """

    @pytest.mark.parametrize(
        "existing, input_line",
        [
            ("vacf.dat", ""),
            ("vacf.log", "log_file = vacf.log\n"),
            ("spectrum.dat", ""),
            ("windowed.dat", "windowed_out_file = windowed.dat\n"),
        ],
    )
    def test_fails_before_trajectory_is_read(
        self, tmp_path, monkeypatch, existing, input_line
    ):
        """
        An existing output, log, spectrum or windowed file raises before
        the trajectory is opened, so no time is lost on large inputs.
        """
        shutil.copytree(EXAMPLES_WATER, tmp_path, dirs_exist_ok=True)
        monkeypatch.chdir(tmp_path)
        with open("vacf.in", "a", encoding="utf-8") as file:
            file.write(input_line)
        (tmp_path / existing).write_text("already there\n", encoding="utf-8")

        vacf_api = importlib.import_module("PQAnalysis.analysis.vacf.api")
        monkeypatch.setattr(vacf_api, "TrajectoryReader", _fail)
        monkeypatch.setattr(vacf_api, "read_restart_file", _fail)

        with pytest.raises(FileWritingModeError, match="already exists"):
            vacf("vacf.in")

        content = (tmp_path / existing).read_text(encoding="utf-8")
        assert content == "already there\n"
