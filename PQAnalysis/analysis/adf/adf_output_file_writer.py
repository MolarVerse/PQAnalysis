"""
A module containing the classes for writing related to an
:py:class:`~PQAnalysis.analysis.adf.adf.ADF` analysis to a file.
"""

# 3rd party imports
from beartype.typing import Tuple

# local imports
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.io import BaseWriter
from PQAnalysis.utils import __header__
from PQAnalysis.type_checking import runtime_type_checking

from .adf import ADF



class ADFDataWriter(BaseWriter):

    """
    Class for writing the data of an
    :py:class:`~PQAnalysis.analysis.adf.adf.ADF`
    analysis to a file.

    The written columns are the bin-center angle (degrees), the
    normalized ADF (a probability density whose integral over the angle
    range is one), the raw angle counts and the sine-corrected ADF (the
    plain ADF divided by ``sin(theta)``, again normalized to unit
    integral).
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
        Np1DNumberArray]
    ):
        """
        Writes the data to the file.

        Parameters
        ----------
        data : Tuple[Np1DNumberArray, Np1DNumberArray,
            Np1DNumberArray, Np1DNumberArray]
            the data output from the ADF.run() method (bin-center angle,
            normalized ADF, raw counts and sine-corrected ADF)
        """
        super().open()

        for i in range(len(data[0])):
            print(
                (
                f"{data[0][i]} {data[1][i]} {data[2][i]} "
                f"{data[3][i]}"
                ),
                file=self.file
            )

        super().close()



class ADFLogWriter(BaseWriter):

    """
    Class for writing the log (setup parameters) of an
    :py:class:`~PQAnalysis.analysis.adf.adf.ADF` analysis
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
    def write_before_run(self, adf: ADF):
        """
        Writes the log before the
        :py:class:`~PQAnalysis.analysis.adf.adf.ADF`
        run() method is called.

        This includes the general header of PQAnalysis and the most
        important setup parameters of the
        :py:class:`~PQAnalysis.analysis.adf.adf.ADF` analysis.

        Parameters
        ----------
        adf : ADF
            the ADF analysis object
        """
        super().open()

        if self.filename is not None:
            print(__header__, file=self.file)
            print(file=self.file)

        print("ADF calculation:", file=self.file)
        print(file=self.file)

        degree = '°'

        # fmt: off
        print(f"    Number of angle bins: {adf.n_angle_bins}", file=self.file)
        print(f"    Angle bin width:      {adf.delta_angle} {degree}", file=self.file)
        print(f"    Angle range:          0 {degree} - {adf.max_angle} {degree}", file=self.file)
        print(file=self.file)
        # fmt: on

        self._write_gates(adf)

        print(f"    Number of frames: {adf.n_frames}", file=self.file)
        print(f"    Number of atoms:  {adf.n_atoms}", file=self.file)
        print(file=self.file)

        print(
            "    Reference selection:",
            adf.reference_selection,
            file=self.file
        )
        print(
            "    total number of atoms in reference selection:",
            len(adf.reference_indices),
            file=self.file
        )
        print(
            "    Target selection (ligand 1):",
            adf.target_selection,
            file=self.file
        )
        print(
            "    total number of atoms in target selection (ligand 1):",
            len(adf.target_indices),
            file=self.file
        )
        print(
            "    Target selection (ligand 2):",
            adf.target_selection_2,
            file=self.file
        )
        print(
            "    total number of atoms in target selection (ligand 2):",
            len(adf.target_indices_2),
            file=self.file
        )
        print(file=self.file)

        print(file=self.file)
        print(file=self.file)
        print(file=self.file)
        print(file=self.file)

        super().close()

    def _write_gates(self, adf: ADF):
        """
        Writes the radial gate settings of the ADF analysis.

        Parameters
        ----------
        adf : ADF
            the ADF analysis object
        """

        angstrom = 'Å'

        if adf.r_min_1 is not None or adf.r_max_1 is not None:
            print(
                f"    i-j radial gate: [{adf.r_min_1}, {adf.r_max_1}) {angstrom}",
                file=self.file
            )
        else:
            print("    i-j radial gate: none", file=self.file)

        if adf.r_min_2 is not None or adf.r_max_2 is not None:
            print(
                f"    i-k radial gate: [{adf.r_min_2}, {adf.r_max_2}) {angstrom}",
                file=self.file
            )
        else:
            print("    i-k radial gate: none", file=self.file)

        print(file=self.file)

    @runtime_type_checking
    def write_after_run(self, adf: ADF):
        """
        Writes the log after the
        :py:class:`~PQAnalysis.analysis.adf.adf.ADF`
        run() method is called.

        This includes the elapsed time of the
        :py:class:`~PQAnalysis.analysis.adf.adf.ADF`
        run() method.

        Parameters
        ----------
        adf : ADF
            the ADF analysis object
        """
        super().open()

        print(f"    Elapsed time: {adf.elapsed_time} s", file=self.file)

        super().close()
