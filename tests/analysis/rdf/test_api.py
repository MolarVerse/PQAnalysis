"""
A module to test the RDF API.
"""

import shutil
from pathlib import Path

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


EXAMPLES_WATER = Path(__file__).resolve().parents[3] / "examples" / "water"


class TestRDFAPIFromOtherDirectory:

    def test_input_file_in_subdirectory(self, tmp_path, monkeypatch):
        """
        The bundled example runs from a parent directory: filenames in the
        input file resolve against the input file, export files against
        the working directory.
        """
        shutil.copytree(EXAMPLES_WATER, tmp_path / "examples" / "water")
        monkeypatch.chdir(tmp_path)

        rdf("examples/water/rdf.in", export_files=["rdf.csv"])

        assert (tmp_path / "examples" / "water" / "rdf.dat").is_file()
        assert (tmp_path / "rdf.csv").is_file()
        assert not (tmp_path / "rdf.dat").exists()
