"""
A module to test the thermal expansion API.
"""

import pytest  # pylint: disable=unused-import

from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.analysis.thermal_expansion.api import thermal_expansion
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception


class TestThermalExpansionAPI:

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
            function=thermal_expansion,
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
            function=thermal_expansion,
            input_file="test",
            md_format=1,
        )

    # @pytest.mark.parametrize("input_file", ["tests/analysis/thermal_expansion/test_input_file.txt"])

    # def test_themal_expansion(self, input_file):
    #     thermal_expansion(input_file=input_file, md_format=MDEngineFormat.PQ)
