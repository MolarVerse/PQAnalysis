"""
A module to test the RDF API.
"""

import pytest  # pylint: disable=unused-import

from PQAnalysis.analysis.rdf.api import rdf
from PQAnalysis.type_checking import _get_type_error_message

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception


class TestRDFAPI:
    def test_wrong_param_types(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=_get_type_error_message(
                "input_file", 1, str,
            ),
            exception=TypeError,
            function=rdf,
            input_file=1,
        )
