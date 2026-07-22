"""
A package containing classes and functions to broaden stick spectra.

Classes
-------
:py:class:`~PQAnalysis.analysis.spectrum_broadening.SpectrumDataWriter`
    A class to write broadened spectra to output files.

Functions
---------
:py:func:`~PQAnalysis.analysis.spectrum_broadening.api.build_spectrum`
    A function to broaden a stick spectrum file and write the result.
:py:func:`~PQAnalysis.analysis.spectrum_broadening.spectrum_broadening.alpha_from_fwhm`
    A function to convert a Gaussian full width at half maximum to
    the exponent alpha.
:py:func:`~PQAnalysis.analysis.spectrum_broadening.spectrum_broadening.fwhm_from_alpha`
    A function to convert a Gaussian exponent alpha to the full width
    at half maximum.
:py:func:`~PQAnalysis.analysis.spectrum_broadening.spectrum_broadening.read_stick_spectrum`
    A function to read a two-column stick spectrum file.
:py:func:`~PQAnalysis.analysis.spectrum_broadening.spectrum_broadening.wavenumber_grid`
    A function to build the regular wavenumber grid of the broadened
    spectrum.
:py:func:`~PQAnalysis.analysis.spectrum_broadening.spectrum_broadening.broaden`
    A function to broaden a stick spectrum on a wavenumber grid.
"""

from .api import build_spectrum
from .spectrum_broadening import (
    DEFAULT_ALPHA,
    DEFAULT_WAVENUMBER_MAX,
    DEFAULT_WAVENUMBER_MIN,
    DEFAULT_WAVENUMBER_STEP,
    alpha_from_fwhm,
    broaden,
    fwhm_from_alpha,
    read_stick_spectrum,
    wavenumber_grid,
)
from .spectrum_broadening_output_file_writer import SpectrumDataWriter

__all__ = [
    "DEFAULT_ALPHA",
    "DEFAULT_WAVENUMBER_MAX",
    "DEFAULT_WAVENUMBER_MIN",
    "DEFAULT_WAVENUMBER_STEP",
    "SpectrumDataWriter",
    "alpha_from_fwhm",
    "broaden",
    "build_spectrum",
    "fwhm_from_alpha",
    "read_stick_spectrum",
    "wavenumber_grid",
]
