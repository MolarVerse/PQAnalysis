import pytest

from PQAnalysis.analysis.bulk_modulus.bulk_modulus_output_file_writer import BulkModulusDataWriter, BulkModulusLogWriter
from PQAnalysis.analysis.bulk_modulus.bulk_modulus import BulkModulus
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

# import topology marker
from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception


class TestBulkModulusLogWriter:

    def test__type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "filename",
                1.0,
                str | None
            ),
            exception=PQTypeError,
            function=BulkModulusLogWriter,
            filename=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("bulk_modulus",
                                                   1.0,
                                                   BulkModulus),
            exception=PQTypeError,
            function=BulkModulusLogWriter("test.out").write_before_run,
            bulk_modulus=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("bulk_modulus",
                                                   1.0,
                                                   BulkModulus),
            exception=PQTypeError,
            function=BulkModulusLogWriter("test.out").write_after_run,
            bulk_modulus=1.0,
            units="GPa",

        )
