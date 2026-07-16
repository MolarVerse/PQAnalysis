"""
A module containing the writer for broadened spectra.
"""

# 3rd party imports
from beartype.typing import Tuple

# local imports
from PQAnalysis.io import BaseWriter
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.type_checking import runtime_type_checking



class SpectrumDataWriter(BaseWriter):

    """
    Class for writing a broadened spectrum to a file.

    Each row contains one grid point and the broadened intensity at
    that grid point in the legacy ``'%8.4f    %16.12e'`` format of
    ``build_spectrum.sh``.
    """

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
    def write(self, data: Tuple[Np1DNumberArray, Np1DNumberArray]) -> None:
        """
        Writes the broadened spectrum to the file.

        Parameters
        ----------
        data : Tuple[Np1DNumberArray, Np1DNumberArray]
            The wavenumber grid and the broadened intensities as
            returned by
            :py:func:`~PQAnalysis.analysis.spectrum_broadening.spectrum_broadening.broaden`.
        """
        super().open()

        for grid_point, intensity in zip(data[0], data[1]):
            print(f"{grid_point:8.4f}    {intensity:16.12e}", file=self.file)

        super().close()
