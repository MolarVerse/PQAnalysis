"""
A module containing the classes for writing related to an
:py:class:`~PQAnalysis.analysis.msd.msd.MSD` analysis to a file.
"""

# 3rd party imports
from beartype.typing import Tuple

# local imports
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.io import BaseWriter
from PQAnalysis.utils import __header__
from PQAnalysis.type_checking import runtime_type_checking

from .msd import MSD



class MSDDataWriter(BaseWriter):

    """
    Class for writing the data of an
    :py:class:`~PQAnalysis.analysis.msd.msd.MSD`
    analysis to a file.

    The output file is written in the legacy Diffcalc format:
    one row per lag index with the columns lag index, MSD in x,
    MSD in y and MSD in z (all in Angstrom^2).
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
            the data output from the MSD.run() method
        """
        super().open()

        lags, msd_x, msd_y, msd_z, _ = data

        for i, lag in enumerate(lags):
            print(
                (
                    f"{int(lag):8d}    {msd_x[i]:12.8f}   "
                    f"{msd_y[i]:12.8f}   {msd_z[i]:12.8f}"
                ),
                file=self.file
            )

        super().close()



class MSDLogWriter(BaseWriter):

    """
    Class for writing the log (setup parameters) of an
    :py:class:`~PQAnalysis.analysis.msd.msd.MSD` analysis
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
    def write_before_run(self, msd: MSD):
        """
        Writes the log before the
        :py:class:`~PQAnalysis.analysis.msd.msd.MSD`
        run() method is called.

        This includes the general header of PQAnalysis
        and the most important setup parameters of the
        :py:class:`~PQAnalysis.analysis.msd.msd.MSD` analysis.

        Parameters
        ----------
        msd : MSD
            the MSD analysis object
        """
        super().open()

        if self.filename is not None:
            print(__header__, file=self.file)
            print(file=self.file)

        print("MSD calculation:", file=self.file)
        print(file=self.file)

        # fmt: off
        print(f"    Window size (frames): {msd.window}", file=self.file)
        print(f"    Origin gap (frames):  {msd.gap}", file=self.file)
        print(f"    Start frame:          {msd.n_start}", file=self.file)
        print(f"    Stop frame:           {msd.stop_frame}", file=self.file)
        print(f"    Number of origins:    {msd.total_origins}", file=self.file)
        print(file=self.file)
        # fmt: on

        print(f"    Number of frames: {msd.n_frames}", file=self.file)
        print(f"    Number of atoms:  {msd.n_atoms}", file=self.file)
        print(file=self.file)

        print(
            "    Target selection:",
            msd.target_selection,
            file=self.file
        )
        print(
            "    total number of atoms in target selection:",
            len(msd.target_indices),
            file=self.file
        )
        print(file=self.file)

        if msd.time_step is not None:
            print(
                f"    Time step:  {msd.time_step} ps",
                file=self.file
            )
            print(
                f"    Fit window: last {msd.fit_window} points",
                file=self.file
            )
            print(file=self.file)

        print(file=self.file)
        print(file=self.file)
        print(file=self.file)
        print(file=self.file)

        super().close()

    @runtime_type_checking
    def write_after_run(self, msd: MSD):
        """
        Writes the log after the
        :py:class:`~PQAnalysis.analysis.msd.msd.MSD`
        run() method is called.

        This includes the elapsed time of the
        :py:class:`~PQAnalysis.analysis.msd.msd.MSD`
        run() method and, if a time step was given, the
        diffusion coefficients obtained from the linear fit
        of the MSD tail.

        Parameters
        ----------
        msd : MSD
            the MSD analysis object
        """
        super().open()

        if msd.fit_results is not None:
            print("    Diffusion coefficients (Einstein relation):",
                file=self.file)
            print(file=self.file)

            for label in ("x", "y", "z", "total"):
                fit = msd.fit_results[label]
                print(
                    (
                        f"    D_{label:<5s} = "
                        f"{fit.diffusion_coefficient:16.8e} +/- "
                        f"{fit.diffusion_coefficient_stderr:16.8e} m^2/s "
                        f"(R^2 = {fit.r_squared:.6f})"
                    ),
                    file=self.file
                )

            print(file=self.file)

        print(f"    Elapsed time: {msd.elapsed_time} s", file=self.file)

        super().close()
