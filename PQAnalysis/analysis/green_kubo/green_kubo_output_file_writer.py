"""
A module containing the classes for writing related to an
:py:class:`~PQAnalysis.analysis.green_kubo.green_kubo.GreenKubo`
analysis to a file.
"""

# 3rd party imports
from beartype.typing import Tuple

# local imports
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.io import BaseWriter
from PQAnalysis.utils import __header__
from PQAnalysis.type_checking import runtime_type_checking

from .green_kubo import GreenKubo



class GreenKuboDataWriter(BaseWriter):

    """
    Class for writing the running-integral data of an
    :py:class:`~PQAnalysis.analysis.green_kubo.green_kubo.GreenKubo`
    analysis to a file.

    The output file contains one row per lag with the columns lag time
    in ps, the un-normalized velocity auto-correlation function Cvv in
    (Angstrom / s)^2 and the running diffusion coefficient D_running in
    m^2 / s.
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
        data: Tuple[Np1DNumberArray, Np1DNumberArray, Np1DNumberArray],
    ):
        """
        Writes the data to the file.

        Parameters
        ----------
        data : Tuple[Np1DNumberArray, Np1DNumberArray, Np1DNumberArray]
            the lag times, the un-normalized velocity auto-correlation
            function and the running diffusion coefficient output from
            the GreenKubo.run() method
        """
        super().open()

        lag_times, cvv, d_running = data

        for time_value, cvv_value, d_value in zip(lag_times, cvv, d_running):
            print(
                f"{time_value:12.6f}    {cvv_value:20.10e}    "
                f"{d_value:20.10e}",
                file=self.file,
            )

        super().close()



class GreenKuboLogWriter(BaseWriter):

    """
    Class for writing the log (setup parameters and plateau result) of
    an
    :py:class:`~PQAnalysis.analysis.green_kubo.green_kubo.GreenKubo`
    analysis to a file.
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
    def write_before_run(self, green_kubo: GreenKubo):
        """
        Writes the log before the
        :py:class:`~PQAnalysis.analysis.green_kubo.green_kubo.GreenKubo`
        run() method is called.

        This includes the general header of PQAnalysis and the most
        important setup parameters of the
        :py:class:`~PQAnalysis.analysis.green_kubo.green_kubo.GreenKubo`
        analysis.

        Parameters
        ----------
        green_kubo : GreenKubo
            the GreenKubo analysis object
        """
        super().open()

        if self.filename is not None:
            print(__header__, file=self.file)
            print(file=self.file)

        print(
            "Green-Kubo self-diffusion coefficient calculation:",
            file=self.file,
        )
        print(file=self.file)

        # fmt: off
        print(f"    Coefficient:          {green_kubo.coefficient}", file=self.file)
        print(f"    Method:               {green_kubo.method}", file=self.file)
        print(f"    Window size (frames): {green_kubo.window_size}", file=self.file)
        print(f"    Origin gap (frames):  {green_kubo.gap}", file=self.file)
        print(f"    Time step:            {green_kubo.time_step} ps", file=self.file)
        print(
            f"    Fit window:           [{green_kubo.fit_start}, "
            f"{green_kubo.fit_stop}] * window",
            file=self.file,
        )
        print(f"    Number of blocks:     {green_kubo.n_blocks}", file=self.file)
        print(file=self.file)
        # fmt: on

        print(f"    Number of frames: {green_kubo.n_frames}", file=self.file)
        print(f"    Number of atoms:  {green_kubo.n_atoms}", file=self.file)
        print(file=self.file)

        print(
            "    Target selection:",
            green_kubo.target_selection,
            file=self.file,
        )
        print(
            "    total number of atoms in target selection:",
            len(green_kubo.target_indices),
            file=self.file,
        )
        print(file=self.file)

        super().close()

    @runtime_type_checking
    def write_after_run(self, green_kubo: GreenKubo):
        """
        Writes the log after the
        :py:class:`~PQAnalysis.analysis.green_kubo.green_kubo.GreenKubo`
        run() method is called.

        This includes the number of time origins, the plateau diffusion
        coefficient with its block-averaged standard error in m^2 / s
        and cm^2 / s, the plateau spread of the running integral and the
        elapsed time of the
        :py:class:`~PQAnalysis.analysis.green_kubo.green_kubo.GreenKubo`
        run() method. The ``+/-`` is the block-averaged standard error
        (over ``n_blocks`` blocks), i.e. an estimate of the run-to-run
        statistical error of ``D``.

        Parameters
        ----------
        green_kubo : GreenKubo
            the GreenKubo analysis object
        """
        super().open()

        print(f"    Number of origins: {green_kubo.n_origins}", file=self.file)
        print(file=self.file)

        n_blocks = len(green_kubo.block_diffusion_coefficients)

        print(
            "    Green-Kubo self-diffusion coefficient (plateau):",
            file=self.file,
        )
        print(file=self.file)
        print(
            (
                f"    D = {green_kubo.diffusion_coefficient:16.8e} +/- "
                f"{green_kubo.diffusion_coefficient_stderr:16.8e} m^2/s "
                f"(n_blocks={n_blocks})"
            ),
            file=self.file,
        )
        print(
            (
                f"    D = "
                f"{green_kubo.diffusion_coefficient_cm2_per_s:16.8e} +/- "
                f"{green_kubo.diffusion_coefficient_stderr_cm2_per_s:16.8e} "
                f"cm^2/s (n_blocks={n_blocks})"
            ),
            file=self.file,
        )
        print(file=self.file)

        print(
            (
                "    Plateau spread (std of running integral over the fit "
                f"window): {green_kubo.diffusion_coefficient_plateau_spread:16.8e}"
                " m^2/s"
            ),
            file=self.file,
        )
        print(file=self.file)

        print(f"    Elapsed time: {green_kubo.elapsed_time} s", file=self.file)

        super().close()
