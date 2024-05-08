"""
A module containing different exceptions and warnings for the :py:class:`~PQAnalysis.analysis.diffcalc.diffcalc.diffcalc` class

"""

from ...exceptions import PQException, PQWarning


class DiffCalcError(PQException):
    """
    Exception raised if something goes wrong during the diffcalc setup or calculation.
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)


class DiffCalcWarning(PQWarning):
    """
    Warning raised if something goes wrong during the diffcalc setup or calculation.
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)