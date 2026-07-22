"""
Tests for the API functions of the VACF analysis.
"""

import os

import numpy as np
import pytest

from PQAnalysis.analysis.vacf import read_static_charges, vacf
from PQAnalysis.analysis.vacf.exceptions import VACFError
from PQAnalysis.io.exceptions import FileWritingModeError

from .. import pytestmark  # pylint: disable=unused-import

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
