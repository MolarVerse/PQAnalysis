"""
A module containing the classes for writing related to an
:py:class:`~PQAnalysis.analysis.finite_differentiation.finite_difference.FiniteDifference` analysis to a file.
"""

# 3rd party imports
from beartype.typing import Tuple

# local imports
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.io import BaseWriter
from PQAnalysis.utils import __header__
from PQAnalysis.type_checking import runtime_type_checking

from .finite_difference import FiniteDifference


class FiniteDifferenceDataWriter(BaseWriter):

    """
    Class for writing the data of an
    :py:class:`~PQAnalysis.analysis.finite_differentiation.finite_difference.FiniteDifference`
    analysis to a file.
    """

    @runtime_type_checking
    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            the filename to write to
        """
        self.filename = filename
        super().__init__(filename)

    @runtime_type_checking
    def write(
        self,
        data: Np1DNumberArray
    ):
        """
        Writes the data to the file.

        Parameters
        ----------
        data : Np1DNumberArray
            the data output from the FiniteDifference.run() method
        """
        super().open()

        for i in range(len(data)):
            self.file.write(f"{data[i]}\n")

        super().close()


class FiniteDifferenceLogWriter(BaseWriter):

    """
    Class for writing the log of an 
    :py:class:`~PQAnalysis.analysis.finite_differentiation.finite_difference.FiniteDifference`
    analysis to a file.
    """

    @runtime_type_checking
    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            the filename to write to
        """
        self.filename = filename
        super().__init__(filename)

    @runtime_type_checking
    def write_before_run(self, finite_difference: FiniteDifference) -> None:
        """
      Writes the log before the 
        :py:class:`~PQAnalysis.analysis.finite_differentiation.finite_difference.FiniteDifference`
        run() method is called.

        This includes the general header of PQAnalysis
        and the most important setup parameters of the
        :py:class:`~PQAnalysis.analysis.finite_differentiation.finite_difference.FiniteDifference` analysis.

        Parameters
        ----------
        finite_difference : FiniteDifference
            the finite difference object
        """
        super().open()

        if self.filename is not None:
            print(__header__, file=self.file)
            print(file=self.file)

        print("Finite difference calculation:", file=self.file)
        print(file=self.file)
        self.file.write(str(finite_difference))

        for i in range(len(finite_difference.finite_difference_points)):
            self.file.write(f"Finite Difference Points: {
                            finite_difference.finite_difference_points[i]}\n")

        super().close()

    @runtime_type_checking
    def write_after_run(self, finite_difference: FiniteDifference) -> None:
        """
        Writes the log after the 
        :py:class:`~PQAnalysis.analysis.finite_differentiation.finite_difference.FiniteDifference` run() method is called.

        This includes the elapsed time of the
        :py:class:`~PQAnalysis.analysis.finite_differentiation.finite_difference.FiniteDifference` run() method.

        Parameters
        ----------
        finite_difference : FiniteDifference
            the finite difference object
        """
        super().open()

        self.file.write(f"Elapsed time: {finite_difference.elapsed_time} s\n")

        super().close()
