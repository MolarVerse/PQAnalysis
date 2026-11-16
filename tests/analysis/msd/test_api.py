"""
Tests for the MSD API.
"""

import importlib
import shutil
from pathlib import Path

import numpy as np
import pytest

from PQAnalysis.analysis.msd.api import msd
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError
from PQAnalysis.io.exceptions import FileWritingModeError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception

EXAMPLES_WATER = Path(__file__).resolve().parents[3] / "examples" / "water"



class TestMSDAPI:

    """
    Tests for the MSD API.
    """

    def test_wrong_param_types(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "input_file",
                1,
                str,
            ),
            exception=PQTypeError,
            function=msd,
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
            function=msd,
            input_file="test",
            md_format=1,
        )

    @pytest.mark.parametrize("example_dir", ["msd"], indirect=False)
    def test_msd(self, test_with_data_dir):
        msd("input.in")

        out_file = Path("msd.dat")
        log_file = Path("msd.log")

        assert out_file.is_file()
        assert log_file.is_file()

        # input.in uses the same window and gap as the legacy
        # reference run diffcalc_O.in
        result = np.loadtxt("msd.dat")
        reference = np.loadtxt("msd_ref_O.dat")

        assert np.array_equal(result[:, 0], reference[:, 0])
        assert np.allclose(
            result[:, 1:],
            reference[:, 1:],
            rtol=1e-5,
            atol=5e-6
        )

        log_contents = log_file.read_text(encoding="utf-8")

        assert "MSD calculation:" in log_contents
        assert "Diffusion coefficients (Einstein relation):" in log_contents
        assert "Elapsed time:" in log_contents



def _fail(*_args, **_kwargs):
    raise AssertionError("trajectory was read")



class TestMSDAPIExistingOutput:

    """
    Tests for existing output files in the msd api function.
    """

    @pytest.mark.parametrize(
        "existing, input_line",
        [
            ("msd.dat", ""),
            ("msd.log", "log_file = msd.log\n"),
        ],
    )
    def test_fails_before_trajectory_is_read(
        self, tmp_path, monkeypatch, existing, input_line
    ):
        """
        An existing output or log file raises before the trajectory is
        opened, so no time is lost on large inputs.
        """
        shutil.copytree(EXAMPLES_WATER, tmp_path, dirs_exist_ok=True)
        monkeypatch.chdir(tmp_path)
        with open("msd.in", "a", encoding="utf-8") as file:
            file.write(input_line)
        (tmp_path / existing).write_text("already there\n", encoding="utf-8")

        msd_api = importlib.import_module("PQAnalysis.analysis.msd.api")
        monkeypatch.setattr(msd_api, "TrajectoryReader", _fail)
        monkeypatch.setattr(msd_api, "read_restart_file", _fail)

        with pytest.raises(FileWritingModeError, match="already exists"):
            msd("msd.in")

        content = (tmp_path / existing).read_text(encoding="utf-8")
        assert content == "already there\n"
