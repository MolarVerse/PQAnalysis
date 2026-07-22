"""
Numerical routines for line broadening of stick spectra.

A stick spectrum is a list of discrete lines, each given by a
wavenumber (in cm^-1) and an intensity. Broadening convolves each
line with a Gaussian (or Lorentzian) kernel on a regular wavenumber
grid using the peak-height convention: the broadened profile of a
single line reaches exactly the line intensity at the line position,
no area normalization is applied. This reproduces the legacy
``build_spectrum.sh`` awk implementation, which calculates
``I(grid) = sum_k I_k * exp(-alpha * (grid - nu_k)^2)`` with a
default ``alpha`` of 0.0025 cm^-2 (FWHM of about 33.3 cm^-1) on a
grid from 10 cm^-1 (inclusive) to 4000 cm^-1 (exclusive) with a step
of 0.25 cm^-1.
"""

from numbers import Real
from pathlib import Path

import numpy as np

from beartype.typing import Tuple

from PQAnalysis.types import Np1DNumberArray, PositiveReal

from .exceptions import SpectrumBroadeningError

#: Default Gaussian exponent alpha in cm^-2 (legacy ``BROAD`` setting).
DEFAULT_ALPHA = 0.0025

#: Default first grid point in cm^-1.
DEFAULT_WAVENUMBER_MIN = 10.0

#: Default (exclusive) last grid point in cm^-1.
DEFAULT_WAVENUMBER_MAX = 4000.0

#: Default grid spacing in cm^-1.
DEFAULT_WAVENUMBER_STEP = 0.25

#: Supported broadening kernels.
KERNELS = ("gaussian", "lorentzian")

#: Maximum number of elements of the outer-difference matrix that are
#: kept in memory at once. Larger problems are processed in grid chunks.
_GRID_CHUNK_SIZE = 4_194_304



def alpha_from_fwhm(fwhm: PositiveReal) -> float:
    """
    Convert a Gaussian full width at half maximum to the exponent alpha.

    The Gaussian kernel is ``exp(-alpha * delta_nu^2)``, so the full
    width at half maximum is ``2 * sqrt(ln(2) / alpha)`` and therefore
    ``alpha = 4 * ln(2) / fwhm^2``.

    Parameters
    ----------
    fwhm : PositiveReal
        The full width at half maximum in cm^-1.

    Returns
    -------
    float
        The Gaussian exponent alpha in cm^-2.

    Raises
    ------
    SpectrumBroadeningError
        If the full width at half maximum is not positive.
    """
    if fwhm <= 0.0:
        raise SpectrumBroadeningError(
            "The full width at half maximum must be positive."
        )

    return 4.0 * np.log(2.0) / float(fwhm)**2



def fwhm_from_alpha(alpha: PositiveReal) -> float:
    """
    Convert a Gaussian exponent alpha to the full width at half maximum.

    Parameters
    ----------
    alpha : PositiveReal
        The Gaussian exponent alpha in cm^-2.

    Returns
    -------
    float
        The full width at half maximum in cm^-1.

    Raises
    ------
    SpectrumBroadeningError
        If alpha is not positive.
    """
    if alpha <= 0.0:
        raise SpectrumBroadeningError("Alpha must be positive.")

    return 2.0 * np.sqrt(np.log(2.0) / float(alpha))



def read_stick_spectrum(
    filename: str
) -> Tuple[Np1DNumberArray, Np1DNumberArray]:
    """
    Read a two-column stick spectrum file.

    The file must contain one line per stick with the wavenumber in
    cm^-1 in the first column and the intensity in the second column.
    Blank lines and lines starting with ``#`` are ignored. Additional
    columns are ignored as well.

    Parameters
    ----------
    filename : str
        The stick spectrum file to read.

    Returns
    -------
    Tuple[Np1DNumberArray, Np1DNumberArray]
        The wavenumbers and intensities of the sticks as float64 arrays.

    Raises
    ------
    SpectrumBroadeningError
        If the file does not exist or contains a malformed line.
    """
    path = Path(filename)
    if not path.is_file():
        raise SpectrumBroadeningError(
            f"Stick spectrum file '{filename}' not found."
        )

    wavenumbers = []
    intensities = []

    with path.open(encoding="utf-8") as file:
        for line_number, line in enumerate(file, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            columns = stripped.split()
            if len(columns) < 2:
                raise SpectrumBroadeningError(
                    f"Line {line_number} of '{filename}' does not "
                    "contain two columns."
                )

            try:
                wavenumbers.append(float(columns[0]))
                intensities.append(float(columns[1]))
            except ValueError as exception:
                raise SpectrumBroadeningError(
                    f"Line {line_number} of '{filename}' contains "
                    "non-numeric data."
                ) from exception

    return (
        np.asarray(wavenumbers, dtype=np.float64),
        np.asarray(intensities, dtype=np.float64),
    )



def wavenumber_grid(
    wavenumber_min: Real = DEFAULT_WAVENUMBER_MIN,
    wavenumber_max: Real = DEFAULT_WAVENUMBER_MAX,
    wavenumber_step: PositiveReal = DEFAULT_WAVENUMBER_STEP,
) -> Np1DNumberArray:
    """
    Build the regular wavenumber grid of the broadened spectrum.

    The grid starts at ``wavenumber_min`` and increases in steps of
    ``wavenumber_step`` while staying strictly below ``wavenumber_max``,
    exactly like the legacy awk loop
    ``for (i = min; i < max; i += step)``. The default grid therefore
    contains 15960 points from 10.0 cm^-1 to 3999.75 cm^-1.

    Parameters
    ----------
    wavenumber_min : Real, optional
        The first grid point in cm^-1, by default 10.0.
    wavenumber_max : Real, optional
        The exclusive upper bound of the grid in cm^-1, by default 4000.0.
    wavenumber_step : PositiveReal, optional
        The grid spacing in cm^-1, by default 0.25.

    Returns
    -------
    Np1DNumberArray
        The wavenumber grid as a float64 array.

    Raises
    ------
    SpectrumBroadeningError
        If the step is not positive or the upper bound is not larger
        than the lower bound.
    """
    if wavenumber_step <= 0.0:
        raise SpectrumBroadeningError(
            "The wavenumber step must be positive."
        )

    if wavenumber_max <= wavenumber_min:
        raise SpectrumBroadeningError(
            "The maximum wavenumber must be larger than the "
            "minimum wavenumber."
        )

    return np.arange(
        float(wavenumber_min),
        float(wavenumber_max),
        float(wavenumber_step),
        dtype=np.float64,
    )



def broaden(
    wavenumbers: Np1DNumberArray,
    intensities: Np1DNumberArray,
    grid: Np1DNumberArray,
    alpha: PositiveReal = DEFAULT_ALPHA,
    kernel: str = "gaussian",
) -> Np1DNumberArray:
    """
    Broaden a stick spectrum on a wavenumber grid.

    Each grid point accumulates the contributions of all sticks using
    the peak-height convention. For the Gaussian kernel the broadened
    spectrum is ``I(g) = sum_k I_k * exp(-alpha * (g - nu_k)^2)``. For
    the Lorentzian kernel the broadened spectrum is
    ``I(g) = sum_k I_k * gamma^2 / ((g - nu_k)^2 + gamma^2)`` where the
    half width at half maximum ``gamma = sqrt(ln(2) / alpha)`` is chosen
    such that both kernels share the same full width at half maximum
    for a given alpha.

    The calculation uses a vectorized outer difference between the grid
    and the stick positions and is chunked along the grid axis to keep
    the memory footprint bounded for large inputs. All accumulation is
    performed in float64.

    Parameters
    ----------
    wavenumbers : Np1DNumberArray
        The stick positions in cm^-1.
    intensities : Np1DNumberArray
        The stick intensities. Must have the same length as the
        stick positions.
    grid : Np1DNumberArray
        The wavenumber grid in cm^-1.
    alpha : PositiveReal, optional
        The Gaussian exponent alpha in cm^-2, by default 0.0025.
    kernel : str, optional
        The broadening kernel, either ``gaussian`` or ``lorentzian``,
        by default ``gaussian``.

    Returns
    -------
    Np1DNumberArray
        The broadened intensities on the grid as a float64 array.

    Raises
    ------
    SpectrumBroadeningError
        If alpha is not positive, the kernel is unknown or the stick
        positions and intensities have different lengths.
    """
    if alpha <= 0.0:
        raise SpectrumBroadeningError("Alpha must be positive.")

    kernel = kernel.lower()
    if kernel not in KERNELS:
        raise SpectrumBroadeningError(
            f"Unknown kernel '{kernel}'. Options are gaussian "
            "and lorentzian."
        )

    wavenumbers = np.asarray(wavenumbers, dtype=np.float64)
    intensities = np.asarray(intensities, dtype=np.float64)
    grid = np.asarray(grid, dtype=np.float64)

    if wavenumbers.shape != intensities.shape:
        raise SpectrumBroadeningError(
            "The stick wavenumbers and intensities must have "
            "the same length."
        )

    result = np.zeros(grid.size, dtype=np.float64)

    if wavenumbers.size == 0:
        return result

    alpha = float(alpha)
    gamma_squared = np.log(2.0) / alpha

    chunk_size = max(1, _GRID_CHUNK_SIZE // wavenumbers.size)

    for start in range(0, grid.size, chunk_size):
        stop = min(start + chunk_size, grid.size)
        delta = grid[start:stop, None] - wavenumbers[None, :]

        if kernel == "gaussian":
            weights = np.exp(-alpha * delta**2)
        else:
            weights = gamma_squared / (delta**2 + gamma_squared)

        result[start:stop] = weights @ intensities

    return result
