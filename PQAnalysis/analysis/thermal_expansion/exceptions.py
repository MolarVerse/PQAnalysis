"""
A module containing different exceptions and warnings for 
the :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion` class
"""

from ...exceptions import PQException, PQWarning


class ThermalExpansionError(PQException):

    """
    Exception raised if something goes wrong during the thermal expansion coefficient setup or calculation.
    """

    def __init__(self, message: str) -> None:
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        self.message = message
        super().__init__(self.message)


class ThermalExpansionWarning(PQWarning):

    """
    Warning raised if something goes wrong during the thermal expansion coefficient setup or calculation.
    """

    def __init__(self, message: str) -> None:
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        self.message = message
        super().__init__(self.message)
