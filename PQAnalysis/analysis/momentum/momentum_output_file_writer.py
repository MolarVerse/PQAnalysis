"""
A module containing the writer for the data of a
:py:class:`~PQAnalysis.analysis.momentum.momentum.Momentum` analysis.
"""

# local imports
from beartype.typing import Sequence

from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.analysis._output_header import format_output_header
from PQAnalysis.analysis.output import (
    AnalysisDataWriter,
    AnalysisTable,
    MOMENTUM_SCHEMA,
)



class MomentumDataWriter(AnalysisDataWriter):

    """
    Class for writing the data of a
    :py:class:`~PQAnalysis.analysis.momentum.momentum.Momentum`
    analysis to a file.

    Each row contains the one-based frame index and the scaled norm
    of the total linear momentum of that frame, reproducing the
    legacy ``equipartition.jl`` output layout.
    """

    schema = MOMENTUM_SCHEMA
    header = format_output_header(schema.title, schema.header_columns)

    @runtime_type_checking
    def __init__(
        self,
        filename: str | None = None,
        mode: str | FileWritingMode = "w",
        export_files: Sequence[str] | None = None,
    ) -> None:
        """
        Parameters
        ----------
        filename : str | None, optional
            The filename to write to. If None, the output is printed
            to stdout, by default None.
        mode : str | FileWritingMode, optional
            The writing mode, by default "w".
        export_files : Sequence[str] | None, optional
            Additional output filenames, by default None.
        """
        super().__init__(filename, mode=mode, export_files=export_files)

    @runtime_type_checking
    def write(self, data: Np1DNumberArray) -> None:
        """
        Writes the momentum norms to the file.

        Parameters
        ----------
        data : Np1DNumberArray
            The scaled momentum norms as returned by the
            :py:meth:`~PQAnalysis.analysis.momentum.momentum.Momentum.run`
            method.
        """
        table = AnalysisTable.from_columns(
            self.schema,
            (range(1, len(data) + 1), data),
        )

        def write_native(file):
            print(self.header, file=file)

            for frame_index, norm in enumerate(data, start=1):
                print(f"{frame_index}  {norm:.12e}", file=file)

        self.write_table(table, write_native)
