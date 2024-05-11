import pytest

from PQAnalysis.io import write, write_box
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError, PQNotImplementedError
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.traj import Trajectory

from . import pytestmark
from ..conftest import assert_logging_with_exception



class TestWriteAPI:

    def test_write_type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_level="ERROR",
            logging_name="TypeChecking",
            message_to_test=get_type_error_message(
            "filename",
            1,
            "str | None"
            ),
            exception=PQTypeError,
            function=write,
            object_to_write=None,
            filename=1,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_level="ERROR",
            logging_name="TypeChecking",
            message_to_test=get_type_error_message(
            "mode",
            1,
            FileWritingMode | str
            ),
            exception=PQTypeError,
            function=write,
            object_to_write=None,
            filename="a",
            mode=1
        )

    def test_write_not_implemented(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_level="ERROR",
            logging_name="io.write_api",
            message_to_test=(
            "Writing object of type <class 'int'> is not implemented yet."
            ),
            exception=PQNotImplementedError,
            function=write,
            object_to_write=1
        )

    def test_write_box_type_checking(self, caplog):

        assert_logging_with_exception(
            caplog=caplog,
            logging_level="ERROR",
            logging_name="TypeChecking",
            message_to_test=get_type_error_message("traj",
            None,
            Trajectory),
            exception=PQTypeError,
            function=write_box,
            traj=None
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_level="ERROR",
            logging_name="TypeChecking",
            message_to_test=get_type_error_message(
            "filename",
            1,
            "str | None"
            ),
            exception=PQTypeError,
            function=write_box,
            traj=Trajectory(),
            filename=1
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_level="ERROR",
            logging_name="TypeChecking",
            message_to_test=get_type_error_message(
            "output_format",
            1,
            "str | None"
            ),
            exception=PQTypeError,
            function=write_box,
            traj=Trajectory(),
            filename="a",
            output_format=1
        )
