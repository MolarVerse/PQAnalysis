"""
Tests for the MSD API.
"""

from pathlib import Path

import numpy as np
import pytest

from PQAnalysis.analysis.msd.api import msd
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception



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
