"""
API functions for spectrum broadening.
"""

from numbers import Real

from beartype.typing import Tuple

from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.types import Np1DNumberArray, PositiveReal

from .exceptions import SpectrumBroadeningError
from .spectrum_broadening import (
    DEFAULT_ALPHA,
    DEFAULT_WAVENUMBER_MAX,
    DEFAULT_WAVENUMBER_MIN,
    DEFAULT_WAVENUMBER_STEP,
    alpha_from_fwhm,
    broaden,
    read_stick_spectrum,
    wavenumber_grid,
)
from .spectrum_broadening_output_file_writer import SpectrumDataWriter



@runtime_type_checking
def build_spectrum(
    input_file: str,
    output: str | None = None,
    alpha: PositiveReal | None = None,
    fwhm: PositiveReal | None = None,
    wavenumber_min: Real = DEFAULT_WAVENUMBER_MIN,
    wavenumber_max: Real = DEFAULT_WAVENUMBER_MAX,
    wavenumber_step: PositiveReal = DEFAULT_WAVENUMBER_STEP,
    kernel: str = "gaussian",
    mode: str | FileWritingMode = "w",
) -> Tuple[Np1DNumberArray, Np1DNumberArray]:
    """
    Broaden a two-column stick spectrum file and write the result.

    Reads a stick spectrum (wavenumber in cm^-1, intensity), broadens
    it on a regular wavenumber grid using the peak-height convention
    and writes one ``'%8.4f    %16.12e'`` row per grid point. This is
    a port of the legacy ``build_spectrum.sh`` awk implementation.

    Parameters
    ----------
    input_file : str
        The two-column stick spectrum file to read.
    output : str | None, optional
        The output file. If None, the output is printed to stdout,
        by default None.
    alpha : PositiveReal | None, optional
        The Gaussian exponent alpha in cm^-2. Mutually exclusive with
        fwhm, by default None, which corresponds to 0.0025 cm^-2 if
        fwhm is not given either.
    fwhm : PositiveReal | None, optional
        The full width at half maximum in cm^-1 as an alternative way
        to specify the broadening width. Mutually exclusive with
        alpha, by default None.
    wavenumber_min : Real, optional
        The first grid point in cm^-1, by default 10.0.
    wavenumber_max : Real, optional
        The exclusive upper bound of the grid in cm^-1,
        by default 4000.0.
    wavenumber_step : PositiveReal, optional
        The grid spacing in cm^-1, by default 0.25.
    kernel : str, optional
        The broadening kernel, either ``gaussian`` or ``lorentzian``,
        by default ``gaussian``.
    mode : str | FileWritingMode, optional
        The writing mode of the output file, by default "w".

    Returns
    -------
    Tuple[Np1DNumberArray, Np1DNumberArray]
        The wavenumber grid and the broadened intensities.

    Raises
    ------
    SpectrumBroadeningError
        If both alpha and fwhm are specified.
    """
    if alpha is not None and fwhm is not None:
        raise SpectrumBroadeningError(
            "The parameters alpha and fwhm are mutually exclusive. "
            "Please specify only one of them."
        )

    if alpha is None:
        alpha = alpha_from_fwhm(fwhm) if fwhm is not None else DEFAULT_ALPHA

    wavenumbers, intensities = read_stick_spectrum(input_file)

    grid = wavenumber_grid(
        wavenumber_min=wavenumber_min,
        wavenumber_max=wavenumber_max,
        wavenumber_step=wavenumber_step,
    )

    broadened = broaden(
        wavenumbers,
        intensities,
        grid,
        alpha=alpha,
        kernel=kernel,
    )

    SpectrumDataWriter(output, mode=mode).write((grid, broadened))

    return grid, broadened
