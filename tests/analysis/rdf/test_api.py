"""
A module to test the RDF API.
"""

import pytest  # pylint: disable=unused-import

from PQAnalysis.analysis.rdf.api import rdf
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception



class TestRDFAPI:

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
            function=rdf,
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
            function=rdf,
            input_file="test",
            md_format=1,
        )
