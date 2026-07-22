"""
A package containing classes and functions to handle velocity and
charge-flux auto-correlation function (VACF) analyses.

Classes
-------
:py:class:`~PQAnalysis.analysis.vacf.vacf.VACF`
    A class to calculate velocity auto-correlation functions.
:py:class:`~PQAnalysis.analysis.vacf.vacf_input_file_reader.VACFInputFileReader`
    A class to read VACF setups from input files.
:py:class:`~PQAnalysis.analysis.vacf.vacf_output_file_writer.VACFDataWriter`
    A class to write VACF data to output files.
:py:class:`~PQAnalysis.analysis.vacf.vacf_output_file_writer.VACFSpectrumDataWriter`
    A class to write VACF spectra to output files.
:py:class:`~PQAnalysis.analysis.vacf.vacf_output_file_writer.VACFWindowedDataWriter`
    A class to write windowed VACF data to output files.
:py:class:`~PQAnalysis.analysis.vacf.vacf_output_file_writer.VACFLogWriter`
    A class to write log files.

Functions
---------
:py:func:`~PQAnalysis.analysis.vacf.api.vacf`
    A function to calculate VACFs from an input file.
:py:func:`~PQAnalysis.analysis.vacf.api.read_static_charges`
    A function to read legacy static charge files.
:py:func:`~PQAnalysis.analysis.vacf.spectrum.vacf_spectrum`
    A function to calculate the legacy cosine-transform spectrum of a
    correlation function.
:py:func:`~PQAnalysis.analysis.vacf.spectrum.apodization_window`
    A function to calculate the legacy apodization windows.
"""

from .api import read_static_charges, vacf
from .vacf import VACF
from .spectrum import WINDOW_FUNCTIONS, apodization_window, vacf_spectrum
from .vacf_input_file_reader import VACFInputFileReader
from .vacf_output_file_writer import (
    VACFDataWriter,
    VACFLogWriter,
    VACFSpectrumDataWriter,
    VACFWindowedDataWriter,
)
from .exceptions import VACFError, VACFWarning

__all__ = [
    "VACF",
    "VACFDataWriter",
    "VACFError",
    "VACFInputFileReader",
    "VACFLogWriter",
    "VACFSpectrumDataWriter",
    "VACFWarning",
    "VACFWindowedDataWriter",
    "WINDOW_FUNCTIONS",
    "apodization_window",
    "read_static_charges",
    "vacf",
    "vacf_spectrum",
]
