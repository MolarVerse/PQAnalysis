import pytest

from PQAnalysis.analysis.rdf.api import rdf

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging


class TestRDFAPI:
    def test_wrong_param_types(self, caplog):
        assert_logging(
            caplog=caplog,
            logging_name=rdf.__name__,
            logging_level="ERROR",
            message_to_test="Input file must be a string",
            function=rdf,
            input_file=1,
        )
