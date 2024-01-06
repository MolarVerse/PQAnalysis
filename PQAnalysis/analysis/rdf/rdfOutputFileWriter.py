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

from .rdf import RadialDistributionFunction
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

    def __init__(self, filename: str,
                 data: Tuple[Np1DNumberArray, Np1DNumberArray,
                             Np1DNumberArray, Np1DNumberArray, Np1DNumberArray]
                 ) -> None:
        """
        It sets the filename and the data to write.

        Parameters
        ----------
        filename : str
            the filename to write to
        data : Tuple[Np1DNumberArray, Np1DNumberArray, Np1DNumberArray, Np1DNumberArray, Np1DNumberArray]
            the data output from the RadialDistributionFunction.run() method
        """
        self.filename = filename
        self.data = data
        super().__init__(filename)

    def write(self):
        """
        Writes the data to the file.
        """
        super().open()

        for i in range(len(self.data[0])):
            print(f"{self.data[0][i]} {self.data[1][i]} {self.data[2][i]} {self.data[3][i]} {self.data[4][i]}",
                  file=self.file)

        super().close()


class RDFLogWriter(BaseWriter):
    """
    Class for writing the log (setup parameters) of an RDF analysis to a file.

    Examples
    --------
    >>> # to write the log before the RDF run() function is called
    >>> RDFLogWriter("rdf.log", rdf).write_before_run()
    >>>
    >>> # to write the log after the RDF run() function is called
    >>> RDFLogWriter("rdf.log", rdf).write_after_run()
    """

    def __inti__(self, filename: str | None, rdf: RadialDistributionFunction) -> None:
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

    def write_before_run(self):
        """
        Writes the log before the RDF run() function is called.

        This includes the general header of PQAnalysis
        and the most important setup parameters of the RDF analysis.
        """
        super().open()

        print(header, file=self.file)
        print(file=self.file)

        print("RDF calculation:", file=self.file)
        print(file=self.file)

        angstrom = u'\u212B'.encode('utf-8')

        # fmt: off
        print(f"    Number of bins: {self.rdf.n_bins}", file=self.file)
        print(f"    Bin width:      {self.rdf.delta_r} {angstrom}", file=self.file)
        print(f"    Minimum radius: {self.rdf.r_min} {angstrom}", file=self.file)
        print(f"    Maximum radius: {self.rdf.r_max} {angstrom}", file=self.file)
        print(file=self.file)
        # fmt: on

        print(f"    Number of frames: {self.rdf.n_frames}", file=self.file)
        print(f"    Number of atoms:  {self.rdf.n_atoms}", file=self.file)
        print(file=self.file)

        # fmt: off
        print(f"    Reference selection: {self.rdf.reference_selection}", file=self.file)
        print(f"    total number of atoms in reference selection: {len(self.rdf.reference_indices)}", file=self.file)
        print(f"    Target selection:    {self.rdf.target_selection}", file=self.file)
        print(f"    total number of atoms in target selection:    {len(self.rdf.target_indices)}", file=self.file)
        # fmt: on

        super().close()

    def write_after_run(self):
        """
        Writes the log after the RDF run() function is called.

        This includes the elapsed time of the RDF run() function.
        """
        super().open()

        print(f"    Elapsed time: {self.rdf.elapsed_time} ms", file=self.file)

        super().close()
