"""
Tests for the exceptions and warnings of the VACF analysis.
"""

from PQAnalysis.analysis.vacf.exceptions import VACFError, VACFWarning
from PQAnalysis.exceptions import PQException, PQWarning

from .. import pytestmark  # pylint: disable=unused-import



class TestVACFExceptions:

    """
    Tests for the VACF exception and warning classes.
    """

    def test_vacf_error_message(self):
        """
        The VACFError stores the message and is a PQException.
        """
        error = VACFError("something went wrong")

        assert error.message == "something went wrong"
        assert str(error) == "something went wrong"
        assert isinstance(error, PQException)

    def test_vacf_warning_message(self):
        """
        The VACFWarning stores the message and is a PQWarning.
        """
        warning = VACFWarning("watch out")

        assert warning.message == "watch out"
        assert str(warning) == "watch out"
        assert isinstance(warning, PQWarning)
        assert isinstance(warning, Warning)
