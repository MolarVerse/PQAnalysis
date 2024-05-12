import pytest

from PQAnalysis.analysis.rdf import RDFDataWriter, RDFLogWriter, RDF
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

# import topology marker
from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception



class TestRDFLogWriter:

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
            function=RDFLogWriter,
            filename=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("rdf",
            1.0,
            RDF),
            exception=PQTypeError,
            function=RDFLogWriter("test.out").write_before_run,
            rdf=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("rdf",
            1.0,
            RDF),
            exception=PQTypeError,
            function=RDFLogWriter("test.out").write_after_run,
            rdf=1.0,
        )
