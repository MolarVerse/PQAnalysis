"""
A module to test the ADF API.
"""

import pytest  # pylint: disable=unused-import

from PQAnalysis.analysis.adf.api import adf
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception



class TestADFAPI:

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
            function=adf,
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
            function=adf,
            input_file="test",
            md_format=1,
        )

    @pytest.mark.parametrize("example_dir", ["adf"], indirect=False)
    def test_end_to_end(self, test_with_data_dir):
        # runs the full wrapper (reader -> ADF -> writers) and checks
        # that the output file has one row per angle bin
        adf("input.in")

        with open("adf.out", "r", encoding="utf-8") as file:
            lines = [line for line in file if line.strip()]

        assert len(lines) == 180
        # four columns: angle, normalized, counts, sin-corrected
        assert all(len(line.split()) == 4 for line in lines)

    @pytest.mark.parametrize("example_dir", ["adf"], indirect=False)
    def test_end_to_end_gated(self, test_with_data_dir):
        adf("input_gated.in")

        with open("adf.out", "r", encoding="utf-8") as file:
            lines = [line for line in file if line.strip()]

        # delta_angle = 2.0 -> 90 bins
        assert len(lines) == 90
