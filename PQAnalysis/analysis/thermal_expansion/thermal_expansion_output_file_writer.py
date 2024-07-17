"""
A module containing the classes for writing related to an
:py:class:`~PQAnalysis.analysis.finite_differentiation.thermal_expansion.ThermalExpansion` analysis to a file.
"""

# 3rd party imports
from beartype.typing import Tuple

# local imports
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.io import BaseWriter
from PQAnalysis.utils import __header__
from PQAnalysis.type_checking import runtime_type_checking

from .thermal_expansion import ThermalExpansion


class ThermalExpansionDataWriter(BaseWriter):

    """
    Class for writing the data of an
    :py:class:`~PQAnalysis.analysis.finite_differentiation.thermal_expansion.ThermalExpansion`
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
        box_av_data: Np1DNumberArray,
        box_std_data: Np1DNumberArray,
        thermal_expansion_data: Np1DNumberArray,

    ):
        """
        Writes the data to the file.

        Parameters
        ----------
        data : Np1DNumberArray
            the data output from the ThermalExpansion.run() method
        """
        super().open()

        for i in range(len(box_av_data)):
            self.file.write(f"{box_av_data[i][0]} {box_av_data[i][1]} {box_av_data[i][2]} {box_std_data[i][0]} {box_std_data[i][1]} {box_std_data[i][2]} {
                            box_std_data[i][3]} {thermal_expansion_data[i][0]} {thermal_expansion_data[i][1]} {thermal_expansion_data[i][2]}\n")

        super().close()


class ThermalExpansionLogWriter(BaseWriter):

    """
    Class for writing the log of an 
    :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion`
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
    def write_before_run(self, thermal_expansion: ThermalExpansion) -> None:
        """
      Writes the log before the 
        :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion`
        run() method is called.

        This includes the general header of PQAnalysis
        and the most important setup parameters of the
        :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion` analysis.

        Parameters
        ----------
        thermal_expansion : ThermalExpansion
            the linear thermal expansion coefficient object
        """
        super().open()

        if self.filename is not None:
            print(__header__, file=self.file)
            print(file=self.file)

        print("Thermal expansion coefficient calculation:", file=self.file)
        print(file=self.file)
        self.file.write(str(thermal_expansion))
        self.file.write(
            "T a_av a_std b_av b_std c_av c_std volume volume_std thermal_expansion_a thermal_expansion_b thermal_expansion_c\n")
        super().close()

    @runtime_type_checking
    def write_after_run(self, thermal_expansion: ThermalExpansion) -> None:
        """
        Writes the log after the 
        :py:class:`~PQAnalysis.analysis.finite_differentiation.thermal_expansion.ThermalExpansion` run() method is called.

        This includes the elapsed time of the
        :py:class:`~PQAnalysis.analysis.finite_differentiation.thermal_expansion.ThermalExpansion` run() method.

        Parameters
        ----------
        thermal_expansion : ThermalExpansion
            the linear thermal expansion coefficient object
        """
        super().open()

        self.file.write(f"Elapsed time: {thermal_expansion.elapsed_time} s\n")

        super().close()
