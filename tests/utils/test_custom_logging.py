"""
Tests for the custom logger.
"""

import logging

import pytest

from PQAnalysis import __package_name__
from PQAnalysis.exceptions import PQException
from PQAnalysis.utils.custom_logging import setup_logger

from . import pytestmark  # pylint: disable=unused-import



class _TestError(PQException):
    """
    An exception type used by the tests below.
    """



class TestCustomLogger:
    """
    Tests that a raised error does not depend on the logging level.
    """

    @pytest.mark.parametrize(
        "level",
        [logging.DEBUG, logging.INFO, logging.ERROR, logging.CRITICAL],
    )
    def test_error_raises_at_every_level(self, level):
        logger = setup_logger(
            logging.getLogger(__package_name__).getChild("ErrorLevelTest")
        )
        logger.setLevel(level)

        with pytest.raises(_TestError):
            logger.error("raised", exception=_TestError)

    def test_error_raises_when_logging_is_silenced(self):
        """
        A logger silenced above CRITICAL still raises, so that raising the
        level cannot turn an error into a silent no-op.
        """
        logger = setup_logger(
            logging.getLogger(__package_name__).getChild("SilencedTest")
        )
        logger.setLevel(logging.CRITICAL + 1)

        with pytest.raises(_TestError):
            logger.error("raised", exception=_TestError)

        with pytest.raises(_TestError):
            logger.critical("raised", exception=_TestError)
