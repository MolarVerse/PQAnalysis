"""
A module containing the classes for writing related to an
:py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` analysis to a file.
"""

# 3rd party imports
from beartype.typing import Tuple

# local imports
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.io import BaseWriter
from PQAnalysis.utils import __header__
from PQAnalysis.type_checking import runtime_type_checking

from .rdf import RDF



class RDFDataWriter(BaseWriter):

    """
    Class for writing the data of an 
    :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF`
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
        data: Tuple[Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray]
    ):
        """
        Writes the data to the file.

        Parameters
        ----------
        data : Tuple[Np1DNumberArray, Np1DNumberArray,
            Np1DNumberArray, Np1DNumberArray, Np1DNumberArray]
            the data output from the RadialDistributionFunction.run() method
        """
        super().open()

        for i in range(len(data[0])):
            print(
                (
                f"{data[0][i]} {data[1][i]} {data[2][i]} "
                f"{data[3][i]} {data[4][i]}"
                ),
                file=self.file
            )

        super().close()



class RDFLogWriter(BaseWriter):

    """
    Class for writing the log (setup parameters) of an 
    :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` analysis
    to a file.
    """

    @runtime_type_checking
    def __init__(self, filename: str | None) -> None:
        """
        Parameters
        ----------
        filename : str | None
            the filename to write to if None, the output is printed to stdout
        """
        self.filename = filename
        super().__init__(filename)

    @runtime_type_checking
    def write_before_run(self, rdf: RDF):
        """
        Writes the log before the 
        :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF`
        run() method is called.

        This includes the general header of PQAnalysis
        and the most important setup parameters of the
        :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` analysis.

        Parameters
        ----------
        rdf : RadialDistributionFunction
            the RDF analysis object
        """
        super().open()

        if self.filename is not None:
            print(__header__, file=self.file)
            print(file=self.file)

        print("RDF calculation:", file=self.file)
        print(file=self.file)

        angstrom = '\u212B'.encode('utf-8')

        # fmt: off
        print(f"    Number of bins: {rdf.n_bins}", file=self.file)
        print(f"    Bin width:      {rdf.delta_r} {angstrom}", file=self.file)
        print(f"    Minimum radius: {rdf.r_min} {angstrom}", file=self.file)
        print(f"    Maximum radius: {rdf.r_max} {angstrom}", file=self.file)
        print(file=self.file)
        # fmt: on

        print(f"    Number of frames: {rdf.n_frames}", file=self.file)
        print(f"    Number of atoms:  {rdf.n_atoms}", file=self.file)
        print(file=self.file)

        print(
            "    Reference selection:",
            rdf.reference_selection,
            file=self.file
        )
        print(
            "    total number of atoms in reference selection:",
            len(rdf.reference_indices),
            file=self.file
        )
        print(
            "    Target selection:   ",
            {rdf.target_selection},
            file=self.file
        )
        print(
            "    total number of atoms in target selection:   ",
            len(rdf.target_indices),
            file=self.file
        )
        print(file=self.file)

        print(
            "    Eliminate intra molecular contributions:",
            rdf.no_intra_molecular,
            file=self.file
        )
        print(file=self.file)

        print(file=self.file)
        print(file=self.file)
        print(file=self.file)
        print(file=self.file)

        super().close()

    @runtime_type_checking
    def write_after_run(self, rdf: RDF):
        """
        Writes the log after the 
        :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF`
        run() method is called.

        This includes the elapsed time of the
        :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF`
        run() method.

        Parameters
        ----------
        rdf : RadialDistributionFunction
            the RDF analysis object
        """
        super().open()

        print(f"    Elapsed time: {rdf.elapsed_time} s", file=self.file)

        super().close()
