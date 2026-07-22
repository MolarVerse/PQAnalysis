"""
Tests for the ADF exceptions and warnings.
"""

from PQAnalysis.analysis.adf.exceptions import ADFError, ADFWarning
from PQAnalysis.exceptions import PQException, PQWarning

from .. import pytestmark  # pylint: disable=unused-import



class TestADFExceptions:

    """
    Tests for the ADF exception and warning classes.
    """

    def test_adf_error(self):
        """The ADF error stores its message and is a PQException."""
        error = ADFError("something went wrong")

        assert error.message == "something went wrong"
        assert str(error) == "something went wrong"
        assert isinstance(error, PQException)

    def test_adf_warning(self):
        """The ADF warning stores its message and is a PQWarning."""
        warning = ADFWarning("watch out")

        assert warning.message == "watch out"
        assert str(warning) == "watch out"
        assert isinstance(warning, PQWarning)
