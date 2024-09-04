import pytest

from PQAnalysis.analysis.thermal_expansion.thermal_expansion_output_file_writer import ThermalExpansionLogWriter, ThermalExpansion
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

# import topology marker
from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception


class TestThermalExpansionLogWriter:

    def test__type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "filename",
                1.0,
                str
            ),
            exception=PQTypeError,
            function=ThermalExpansionLogWriter,
            filename=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("thermal_expansion",
                                                   1.0,
                                                   ThermalExpansion),
            exception=PQTypeError,
            function=ThermalExpansionLogWriter("test.out").write_before_run,
            thermal_expansion=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("thermal_expansion",
                                                   1.0,
                                                   ThermalExpansion),
            exception=PQTypeError,
            function=ThermalExpansionLogWriter("test.out").write_after_run,
            thermal_expansion=1.0,
        )
