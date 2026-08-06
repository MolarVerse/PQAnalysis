"""
A module containing the writer for the data of a
:py:class:`~PQAnalysis.analysis.momentum.momentum.Momentum` analysis.
"""

# local imports
from PQAnalysis.io import BaseWriter
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.analysis._output_header import format_output_header



class MomentumDataWriter(BaseWriter):

    """
    Class for writing the data of a
    :py:class:`~PQAnalysis.analysis.momentum.momentum.Momentum`
    analysis to a file.

    Each row contains the one-based frame index and the scaled norm
    of the total linear momentum of that frame, reproducing the
    legacy ``equipartition.jl`` output layout.
    """

    header = format_output_header(
        "Total linear momentum",
        (
            ("n", "1"),
            ("s*|P(n)|", "scale-dependent"),
        ),
    )

    @runtime_type_checking
    def __init__(
        self,
        filename: str | None = None,
        mode: str | FileWritingMode = "w",
    ) -> None:
        """
        Parameters
        ----------
        filename : str | None, optional
            The filename to write to. If None, the output is printed
            to stdout, by default None.
        mode : str | FileWritingMode, optional
            The writing mode, by default "w".
        """
        self.filename = filename
        super().__init__(filename, mode=mode)

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
        super().open()

        print(self.header, file=self.file)

        for frame_index, norm in enumerate(data, start=1):
            print(f"{frame_index}  {norm:.12e}", file=self.file)

        super().close()
