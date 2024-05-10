import pytest
import os
import shutil
import logging

from contextlib import contextmanager
from _pytest.logging import LogCaptureHandler, _remove_ansi_escape_sequences
from beartype.roar import BeartypeCallHintParamViolation

from PQAnalysis.utils.custom_logging import CustomLoggerException
from PQAnalysis import __package_name__

from . import __beartype_level__


@pytest.fixture(scope="function")
def tmpdir():

    tmpdir = "tmpdir"

    if os.path.exists(tmpdir) and os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    os.mkdir(tmpdir)

    os.chdir(tmpdir)

    yield tmpdir

    os.chdir("..")
    shutil.rmtree(tmpdir)


@pytest.fixture(scope="function")
def test_with_data_dir(example_dir):

    tmpdir = "tmpdir"

    if os.path.exists(tmpdir) and os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)

    shutil.copytree(os.path.join("tests/data", example_dir), tmpdir)

    os.chdir(tmpdir)

    yield tmpdir

    os.chdir("..")
    shutil.rmtree(tmpdir)


class CatchLogFixture:
    """
    Fixture to capture logs regardless of the Propagate flag. See
    https://github.com/pytest-dev/pytest/issues/3697 for details.
    """

    @property
    def text(self) -> str:
        return _remove_ansi_escape_sequences(self.handler.stream.getvalue())

    @contextmanager
    def catch_logs(self, level: int, logger: logging.Logger) -> LogCaptureHandler:
        """Set the level for capturing of logs. After the end of the 'with' statement,
        the level is restored to its original value.
        """
        self.handler = LogCaptureHandler()
        orig_level = logger.level
        orig_propagate = logger.propagate
        logger.setLevel(level)
        logger.addHandler(self.handler)
        logger.propagate = True
        try:
            yield self
        finally:
            logger.setLevel(orig_level)
            logger.removeHandler(self.handler)
            logger.propagate = orig_propagate


@pytest.fixture
def capture_log():
    return CatchLogFixture().catch_logs


@contextmanager
def caplog_for_logger(caplog, logger_name, level):
    caplog.clear()
    caplog.set_level(level)
    logger = logging.getLogger(logger_name)
    logger.addHandler(caplog.handler)
    yield
    logger.removeHandler(caplog.handler)


def assert_logging_with_exception(caplog, logging_name, logging_level, message_to_test, exception, function, *args, **kwargs):
    with caplog_for_logger(caplog, __package_name__ + "." + logging_name, logging_level):
        result = None
        try:
            result = function(*args, **kwargs)
        except (CustomLoggerException, BeartypeCallHintParamViolation) as e:
            if isinstance(e, BeartypeCallHintParamViolation) and exception is TypeError:
                return result
            if not isinstance(e, CustomLoggerException):
                raise e

        record = caplog.records[0]

        assert record.name == __package_name__ + \
            "." + logging_name
        assert record.levelno == logging.getLevelName(
            logging_level
        )

        if exception is not None:
            assert record.custom_exception == exception

        message = record.msg

        messages = [message.strip() for message in message.split("\n")]

        for message in messages:
            assert message in message_to_test

        return result


def assert_logging(caplog, logging_name, logging_level, message_to_test, function, *args, **kwargs):
    return assert_logging_with_exception(
        caplog,
        logging_name,
        logging_level,
        message_to_test,
        None,
        function,
        *args,
        **kwargs
    )


def assert_type_error_in_debug_mode(func, *args, **kwargs):
    if __beartype_level__ == "DEBUG":
        with pytest.raises(BeartypeCallHintParamViolation):
            func(*args, **kwargs)
    else:
        try:
            func(*args, **kwargs)
        except:
            pytest.fail("Function raised an exception in non-debug mode")
