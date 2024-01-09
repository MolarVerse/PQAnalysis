"""
A module containing the classes for writing related to an RDF analysis to a file.

Classes
-------
RDFDataWriter
    A class for writing the data of an RDF analysis to a file.
RDFLogWriter
    A class for writing the log (setup parameters) of an RDF analysis to a file.
"""

from beartype.typing import Tuple

from .rdf import RDF
from ...types import Np1DNumberArray
from ...io import BaseWriter
from ...utils import header


class RDFDataWriter(BaseWriter):
    """
    Class for writing the data of an RDF analysis to a file.

    Examples
    --------
    >>> RDFDataWriter("rdf.dat", rdf_data).write()
    """

    def __init__(self, filename: str) -> None:
        """
        It sets the filename and the data to write.

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
            the data output from the RadialDistributionFunction.run() method
        """
        super().open()

        for i in range(len(data[0])):
            print(
                f"{data[0][i]} {data[1][i]} {data[2][i]} {data[3][i]} {data[4][i]}", file=self.file)

        super().close()


class RDFLogWriter(BaseWriter):
    """
    Class for writing the log (setup parameters) of an RDF analysis to a file.
    """

    def __init__(self, filename: str | None) -> None:
        """
        It sets the filename and the RDF analysis object.

        Parameters
        ----------
        filename : str | None
            the filename to write to if None, the output is printed to stdout
        rdf : RadialDistributionFunction
            the RDF analysis object
        """
        self.filename = filename
        super().__init__(filename)

    def write_before_run(self, rdf: RDF):
        """
        Writes the log before the RDF run() function is called.

        This includes the general header of PQAnalysis
        and the most important setup parameters of the RDF analysis.

        Parameters
        ----------
        rdf : RadialDistributionFunction
            the RDF analysis object
        """
        super().open()

        if self.filename is not None:
            print(header, file=self.file)
            print(file=self.file)

        print("RDF calculation:", file=self.file)
        print(file=self.file)

        angstrom = u'\u212B'.encode('utf-8')

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

        # fmt: off
        print(f"    Reference selection: {rdf.reference_selection}", file=self.file)
        print(f"    total number of atoms in reference selection: {len(rdf.reference_indices)}", file=self.file)
        print(f"    Target selection:    {rdf.target_selection}", file=self.file)
        print(f"    total number of atoms in target selection:    {len(rdf.target_indices)}", file=self.file)
        print(file=self.file)
        # fmt: on

        #fmt: off
        print(f"    Eliminate intra molecular contributions: {rdf.no_intra_molecular}", file=self.file)
        print(file=self.file)
        #fmt: on

        print(file=self.file)
        print(file=self.file)
        print(file=self.file)
        print(file=self.file)

        super().close()

    def write_after_run(self, rdf: RDF):
        """
        Writes the log after the RDF run() function is called.

        This includes the elapsed time of the RDF run() function.

        Parameters
        ----------
        rdf : RadialDistributionFunction
            the RDF analysis object
        """
        super().open()

        print(f"    Elapsed time: {rdf.elapsed_time} s", file=self.file)

        super().close()
