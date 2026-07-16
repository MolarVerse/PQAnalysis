"""
A module containing the classes for writing related to an
:py:class:`~PQAnalysis.analysis.vacf.vacf.VACF` analysis to a file.
"""

# 3rd party imports
from beartype.typing import Tuple

# local imports
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.io import BaseWriter
from PQAnalysis.utils import __header__
from PQAnalysis.type_checking import runtime_type_checking

from .vacf import VACF



class VACFDataWriter(BaseWriter):

    """
    Class for writing the data of an
    :py:class:`~PQAnalysis.analysis.vacf.vacf.VACF`
    analysis to a file.

    The output file is written in the legacy FreqCalc/Fluxfreqcalc
    format: one row per lag with the columns lag time in ps and
    normalized VACF.
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
    def write(self, data: Tuple[Np1DNumberArray, Np1DNumberArray]):
        """
        Writes the data to the file.

        Parameters
        ----------
        data : Tuple[Np1DNumberArray, Np1DNumberArray]
            the time axis and correlation function output from the
            VACF.run() method
        """
        super().open()

        time, correlation = data

        for time_value, correlation_value in zip(time, correlation):
            print(
                f"{time_value:10.6f}    {correlation_value:12.8f}",
                file=self.file
            )

        super().close()



class VACFSpectrumDataWriter(BaseWriter):

    """
    Class for writing the spectrum of an
    :py:class:`~PQAnalysis.analysis.vacf.vacf.VACF`
    analysis to a file.

    The output file is written in the legacy ft.f format: one row per
    frequency index with the columns wavenumber in cm^-1 and amplitude.
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
    def write(self, data: Tuple[Np1DNumberArray, Np1DNumberArray]):
        """
        Writes the spectrum to the file.

        Parameters
        ----------
        data : Tuple[Np1DNumberArray, Np1DNumberArray]
            the wavenumbers and amplitudes output from the
            :py:func:`~PQAnalysis.analysis.vacf.spectrum.vacf_spectrum`
            function
        """
        super().open()

        wavenumbers, amplitudes = data

        for wavenumber, amplitude in zip(wavenumbers, amplitudes):
            print(
                f"{wavenumber:13.7f}  {amplitude:14.10f}",
                file=self.file
            )

        super().close()



class VACFWindowedDataWriter(BaseWriter):

    """
    Class for writing the windowed (apodized) correlation function of
    an :py:class:`~PQAnalysis.analysis.vacf.vacf.VACF` analysis to a
    file.

    The output file is written in the legacy ft.f windowfile format:
    one row per lag with the columns lag time in ps and windowed
    correlation function.
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
    def write(self, data: Tuple[Np1DNumberArray, Np1DNumberArray]):
        """
        Writes the windowed correlation function to the file.

        Parameters
        ----------
        data : Tuple[Np1DNumberArray, Np1DNumberArray]
            the time axis and the windowed correlation function output
            from the
            :py:func:`~PQAnalysis.analysis.vacf.spectrum.vacf_spectrum`
            function
        """
        super().open()

        time, correlation = data

        for time_value, correlation_value in zip(time, correlation):
            print(
                f"{time_value:9.4f}  {correlation_value:14.10f}",
                file=self.file
            )

        super().close()



class VACFLogWriter(BaseWriter):

    """
    Class for writing the log (setup parameters) of an
    :py:class:`~PQAnalysis.analysis.vacf.vacf.VACF` analysis
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
    def write_before_run(self, vacf: VACF):
        """
        Writes the log before the
        :py:class:`~PQAnalysis.analysis.vacf.vacf.VACF`
        run() method is called.

        This includes the general header of PQAnalysis
        and the most important setup parameters of the
        :py:class:`~PQAnalysis.analysis.vacf.vacf.VACF` analysis.

        Parameters
        ----------
        vacf : VACF
            the VACF analysis object
        """
        super().open()

        if self.filename is not None:
            print(__header__, file=self.file)
            print(file=self.file)

        if vacf.flux:
            print("Charge-flux auto-correlation calculation:",
                  file=self.file)
        else:
            print("VACF calculation:", file=self.file)
        print(file=self.file)

        # fmt: off
        print(f"    Window size (frames): {vacf.window_size}", file=self.file)
        print(f"    Origin gap (frames):  {vacf.gap}", file=self.file)
        print(f"    Stop frame:           {vacf.stop_frame}", file=self.file)
        print(f"    Time step:            {vacf.time_step} ps", file=self.file)
        print(f"    Method:               {vacf.method}", file=self.file)
        print(file=self.file)
        # fmt: on

        print(f"    Number of frames: {vacf.n_frames}", file=self.file)
        print(f"    Number of atoms:  {vacf.n_atoms}", file=self.file)
        print(file=self.file)

        print(
            "    Target selection:",
            vacf.target_selection,
            file=self.file
        )
        print(
            "    total number of atoms in target selection:",
            len(vacf.target_indices),
            file=self.file
        )
        print(file=self.file)

        super().close()

    @runtime_type_checking
    def write_after_run(self, vacf: VACF):
        """
        Writes the log after the
        :py:class:`~PQAnalysis.analysis.vacf.vacf.VACF`
        run() method is called.

        This includes the number of time origins and the elapsed time
        of the :py:class:`~PQAnalysis.analysis.vacf.vacf.VACF`
        run() method.

        Parameters
        ----------
        vacf : VACF
            the VACF analysis object
        """
        super().open()

        print(f"    Number of origins: {vacf.n_origins}", file=self.file)
        print(file=self.file)

        print(f"    Elapsed time: {vacf.elapsed_time} s", file=self.file)

        super().close()
