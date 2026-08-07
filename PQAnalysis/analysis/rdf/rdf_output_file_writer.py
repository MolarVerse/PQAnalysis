"""
A module containing the classes for writing related to an
:py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` analysis to a file.
"""

# 3rd party imports
from beartype.typing import Sequence, Tuple

# local imports
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.io import BaseWriter
from PQAnalysis.utils import __header__
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.analysis._output_header import format_output_header
from PQAnalysis.analysis.output import (
    AnalysisDataWriter,
    AnalysisTable,
    RDF_SCHEMA,
)

from .rdf import RDF



class RDFDataWriter(AnalysisDataWriter):

    """
    Class for writing the data of an 
    :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF`
    analysis to a file.

    Each row contains five columns: bin-center distance in Angstrom,
    radial distribution function, cumulative coordination number,
    density-normalized shell population in Angstrom^3 and ideal-gas
    pair-count residual. See :ref:`RDF output files <analysis-output-rdf>`
    for the exact definitions and normalization formulas.
    """

    schema = RDF_SCHEMA
    header = format_output_header(schema.title, schema.header_columns)

    @runtime_type_checking
    def __init__(
        self,
        filename: str,
        export_files: Sequence[str] | None = None,
    ) -> None:
        """
        Parameters
        ----------
        filename : str
            The primary output filename. Its extension selects the format.
        export_files : Sequence[str] | None, optional
            Additional output filenames, by default None.
        """
        super().__init__(filename, export_files=export_files)

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
            The bin centers, radial distribution function, cumulative
            coordination number, density-normalized shell population and
            ideal-gas pair-count residual returned by
            :py:meth:`~PQAnalysis.analysis.rdf.rdf.RDF.run`.
        """
        table = AnalysisTable.from_columns(self.schema, data)

        def write_native(file):
            print(self.header, file=file)

            for i in range(len(data[0])):
                print(
                    (
                        f"{data[0][i]} {data[1][i]} {data[2][i]} "
                        f"{data[3][i]} {data[4][i]}"
                    ),
                    file=file
                )

        self.write_table(table, write_native)



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
        print("    Target selection:   ", rdf.target_selection, file=self.file)
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
