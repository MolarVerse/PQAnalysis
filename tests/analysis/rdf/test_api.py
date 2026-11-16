"""
A module to test the RDF API.
"""

import shutil
from pathlib import Path

import pytest  # pylint: disable=unused-import

from PQAnalysis.analysis.rdf import api as rdf_api
from PQAnalysis.analysis.rdf.api import rdf
from PQAnalysis.io.exceptions import FileWritingModeError
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception

EXAMPLES_WATER = Path(__file__).resolve().parents[3] / "examples" / "water"


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


class TestRDFAPIExistingOutput:

    @pytest.mark.parametrize(
        "existing, input_line",
        [
            ("rdf.dat", ""),
            ("rdf.log", "log_file = rdf.log\n"),
        ],
    )
    def test_fails_before_trajectory_is_read(
        self, tmp_path, monkeypatch, existing, input_line
    ):
        """
        An existing output or log file raises before the trajectory is
        opened, so no time is lost on large inputs.
        """
        shutil.copytree(EXAMPLES_WATER, tmp_path, dirs_exist_ok=True)
        monkeypatch.chdir(tmp_path)
        with open("rdf.in", "a", encoding="utf-8") as file:
            file.write(input_line)
        (tmp_path / existing).write_text("already there\n", encoding="utf-8")

        def fail(*_args, **_kwargs):
            raise AssertionError("trajectory was read")

        monkeypatch.setattr(rdf_api, "TrajectoryReader", fail)
        monkeypatch.setattr(rdf_api, "RestartFileReader", fail)

        with pytest.raises(FileWritingModeError, match="already exists"):
            rdf("rdf.in")

        assert (tmp_path / existing).read_text(encoding="utf-8") == (
            "already there\n"
        )
