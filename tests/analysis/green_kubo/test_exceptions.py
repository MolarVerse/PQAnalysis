"""
Tests for the Green-Kubo exceptions and warnings.
"""

from PQAnalysis.analysis.green_kubo.exceptions import (
    GreenKuboError,
    GreenKuboWarning,
)
from PQAnalysis.exceptions import PQException, PQWarning

from .. import pytestmark  # pylint: disable=unused-import



class TestGreenKuboExceptions:

    """
    Tests for the Green-Kubo exception and warning classes.
    """

    def test_green_kubo_error(self):
        error = GreenKuboError("something went wrong")

        assert error.message == "something went wrong"
        assert str(error) == "something went wrong"
        assert isinstance(error, PQException)

    def test_green_kubo_warning(self):
        warning = GreenKuboWarning("watch out")

        assert warning.message == "watch out"
        assert str(warning) == "watch out"
        assert isinstance(warning, PQWarning)
