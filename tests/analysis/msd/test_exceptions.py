"""
Tests for the MSD exceptions and warnings.
"""

from PQAnalysis.analysis.msd.exceptions import MSDError, MSDWarning
from PQAnalysis.exceptions import PQException, PQWarning

from .. import pytestmark  # pylint: disable=unused-import



class TestMSDExceptions:

    """
    Tests for the MSD exception and warning classes.
    """

    def test_msd_error(self):
        error = MSDError("something went wrong")

        assert error.message == "something went wrong"
        assert str(error) == "something went wrong"
        assert isinstance(error, PQException)

    def test_msd_warning(self):
        warning = MSDWarning("watch out")

        assert warning.message == "watch out"
        assert str(warning) == "watch out"
        assert isinstance(warning, PQWarning)
