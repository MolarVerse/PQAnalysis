"""
A module containing the classes for writing related to an :py:class:`~PQAnalysis.analysis.diffcalc.diffcalc.diffcalc` analysis to a file.
"""

# 3rd party imports
from beartype.typing import Tuple

# local imports
from .diffcalc import DiffCalc
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.io import BaseWriter
from PQAnalysis.utils import header


class DiffCalcDataWriter(BaseWriter):
    """
    Class for writing the data of an :py:class:`~PQAnalysis.analysis.diffcalc.diffcalc.diffcalc` analysis to a file.
    """

    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            the filename to write to
        """
        self.filename = filename
        super().__init__(filename)

    def write(self,
              data: Tuple[Np1DNumberArray, Np1DNumberArray,
                          Np1DNumberArray, Np1DNumberArray, Np1DNumberArray]
              ):
        """
        Writes the data to the file.

        Parameters
        ----------
        data : Tuple[Np1DNumberArray, Np1DNumberArray, Np1DNumberArray, Np1DNumberArray, Np1DNumberArray]
            the data output from the SelfDiffusionCoefficientFunction.run() method
        """
        super().open()

        for i in range(len(data[0])):
            print(
                f"{data[0][i]} {data[1][i]} {data[2][i]} {data[3][i]} {data[4][i]}", file=self.file)

        super().close()


class DiffCalcLogWriter(BaseWriter):
    """
    Class for writing the log (setup parameters) of an :py:class:`~PQAnalysis.analysis.diffcalc.diffcalc.diffcalc` analysis to a file.
    """

    def __init__(self, filename: str | None) -> None:
        """
        Parameters
        ----------
        filename : str | None
            the filename to write to if None, the output is printed to stdout
        """
        self.filename = filename
        super().__init__(filename)

    def write_before_run(self, diffcalc: DiffCalc):
        """
        Writes the log before the :py:class:`~PQAnalysis.analysis.diffcalc.diffcalc.diffcalc` run() method is called.

        This includes the general header of PQAnalysis
        and the most important setup parameters of the :py:class:`~PQAnalysis.analysis.diffcalc.diffcalc.diffcalc` analysis.

        Parameters
        ----------
        diffcalc : SelfDiffusionCoefficientFunction
            the diffcalc analysis object
        """
        super().open()

        if self.filename is not None:
            print(header, file=self.file)
            print(file=self.file)

        print("diffcalc calculation:", file=self.file)
        print(file=self.file)

        angstrom = u'\u212B'.encode('utf-8')

        # fmt: off
        print(f"    Number of bins: {diffcalc.n_bins}", file=self.file)
        print(f"    Bin width:      {diffcalc.delta_r} {angstrom}", file=self.file)
        print(f"    Minimum radius: {diffcalc.r_min} {angstrom}", file=self.file)
        print(f"    Maximum radius: {diffcalc.r_max} {angstrom}", file=self.file)
        print(file=self.file)
        # fmt: on

        print(f"    Number of frames: {diffcalc.n_frames}", file=self.file)
        print(f"    Number of atoms:  {diffcalc.n_atoms}", file=self.file)
        print(file=self.file)

        # fmt: off
        print(f"    Reference selection: {diffcalc.reference_selection}", file=self.file)
        print(f"    total number of atoms in reference selection: {len(diffcalc.reference_indices)}", file=self.file)
        print(f"    Target selection:    {diffcalc.target_selection}", file=self.file)
        print(f"    total number of atoms in target selection:    {len(diffcalc.target_indices)}", file=self.file)
        print(file=self.file)
        # fmt: on

        #fmt: off
        print(f"    Eliminate intra molecular contributions: {diffcalc.no_intra_molecular}", file=self.file)
        print(file=self.file)
        #fmt: on

        print(file=self.file)
        print(file=self.file)
        print(file=self.file)
        print(file=self.file)

        super().close()

    def write_after_run(self, diffcalc: DiffCalc):
        """
        Writes the log after the :py:class:`~PQAnalysis.analysis.diffcalc.diffcalc.diffcalc` run() method is called.

        This includes the elapsed time of the :py:class:`~PQAnalysis.analysis.diffcalc.diffcalc.diffcalc` run() method.

        Parameters
        ----------
        diffcalc : SelfDiffusionCoefficientFunction
            the diffcalc analysis object
        """
        super().open()

        print(f"    Elapsed time: {diffcalc.elapsed_time} s", file=self.file)

        super().close()