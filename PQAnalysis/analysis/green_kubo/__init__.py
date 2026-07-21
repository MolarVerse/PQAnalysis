"""
A package containing classes and functions to handle Green-Kubo
transport coefficient analyses. Currently the self-diffusion
coefficient is implemented, obtained from the time integral of the
un-normalized velocity auto-correlation function.

Classes
-------
:py:class:`~PQAnalysis.analysis.green_kubo.green_kubo.GreenKubo`
    A class to calculate Green-Kubo transport coefficients.
:py:class:`~PQAnalysis.analysis.green_kubo.green_kubo_input_file_reader.GreenKuboInputFileReader`
    A class to read Green-Kubo setups from input files.
:py:class:`~PQAnalysis.analysis.green_kubo.green_kubo_output_file_writer.GreenKuboDataWriter`
    A class to write Green-Kubo running-integral data to output files.
:py:class:`~PQAnalysis.analysis.green_kubo.green_kubo_output_file_writer.GreenKuboLogWriter`
    A class to write log files.

Functions
---------
:py:func:`~PQAnalysis.analysis.green_kubo.api.green_kubo`
    A function to calculate Green-Kubo coefficients from an input file.
:py:func:`~PQAnalysis.analysis.green_kubo.green_kubo.velocity_acf_fft`
    A function to calculate the un-normalized velocity auto-correlation
    function with the Wiener-Khinchin (FFT) estimator.
:py:func:`~PQAnalysis.analysis.green_kubo.green_kubo.velocity_acf_direct`
    A function to calculate the un-normalized velocity auto-correlation
    function with a direct sliding-time-origin estimator.
:py:func:`~PQAnalysis.analysis.green_kubo.green_kubo.cumulative_trapezoid`
    A function to calculate a cumulative trapezoidal integral.
"""

from .api import green_kubo
from .green_kubo import (
    GreenKubo,
    ANGSTROM2_PER_S2_PS_TO_M2_PER_S,
    M2_PER_S_TO_CM2_PER_S,
    cumulative_trapezoid,
    velocity_acf_direct,
    velocity_acf_fft,
)
from .green_kubo_input_file_reader import GreenKuboInputFileReader
from .green_kubo_output_file_writer import (
    GreenKuboDataWriter,
    GreenKuboLogWriter,
)
from .exceptions import GreenKuboError, GreenKuboWarning

__all__ = [
    "ANGSTROM2_PER_S2_PS_TO_M2_PER_S",
    "GreenKubo",
    "GreenKuboDataWriter",
    "GreenKuboError",
    "GreenKuboInputFileReader",
    "GreenKuboLogWriter",
    "GreenKuboWarning",
    "M2_PER_S_TO_CM2_PER_S",
    "cumulative_trapezoid",
    "green_kubo",
    "velocity_acf_direct",
    "velocity_acf_fft",
]
