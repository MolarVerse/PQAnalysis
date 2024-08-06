"""
A module containing the classes for finite differentiation analysis.
"""

from beartype.typing import List
import logging

# 3rd party imports
import numpy as np
from scipy.optimize import curve_fit
# 3rd party imports


# local absolute imports
# from PQAnalysis.config import with_progress_bar

from PQAnalysis.physical_data import Energy
from PQAnalysis.types import Np1DNumberArray, Np2DNumberArray, PositiveReal
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__

from PQAnalysis.type_checking import runtime_type_checking
# local relative imports
from .exceptions import BulkModulusError, BulkModulusWarning


def _third_order_bmeos(v, b0, b0_prime, v0):
    """
    Calculate the pressure using the third order Birch-Murnaghan equation of state.

    :math:
        P = 3/2 * B0 * ((V0/V)**(7/3) - (V0/V)**(5/3)) * \
                        (1 + 0.75 * (B0_prime - 4) * ((V0/V)**(2/3) - 1))

    Parameters
    ----------
    v : float
        the volume
    v0 : float
        the equilibrium volume
    b0 : float
        the bulk modulus
    b0_prime : float
        the first derivative of the bulk modulus

    Returns
    -------
    float
        the pressure
    """
    return 3/2 * b0 * ((v0/v)**(7/3) - (v0/v)**(5/3)) * (1 + 0.75 * (b0_prime - 4) * ((v0/v)**(2/3) - 1))


class BulkModulus:
    """
    A class to perform bulk modulus analysis.

    The bulk modulus is calculated using the following formula:

    .. math::
        B = -V \\left(\\frac{dP}{dV}\\right)_T

    .. math::
        \\left(\\frac{dP}{dV}\\right)_T = \\frac{P(V+\\Delta V) - P(V-\\Delta V)}{2\\Delta V}
    -----
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    _mode_default = "simple"

    @runtime_type_checking
    def __init__(
        self,
        pressures_perturbation: Np1DNumberArray | Np2DNumberArray | None = None,
        volumes_perturbation: Np1DNumberArray | None = None,
        volume_equilibrium: PositiveReal | None = None,
        mode: str | None = None,
    ):
        """
        Parameters
        ----------
        pressures_perturbation : Np1DNumberArray | Np2DNumberArray| None
            the averaged and standard deviation (optional)
            pressure data of the perturbation
        volumes_perturbation : Np1DNumberArray | None
            the small perturbation in the volume
        volume_equilibrium : PositiveReal | None
            the equilibrium volume
        mode : str | None
            the mode of the analysis
        """

        if pressures_perturbation is None:
            self.logger.error(
                "The pressure perturbation must be provided.",
                exception=BulkModulusError,
            )
        if volumes_perturbation is None:
            self.logger.error(
                "The volume perturbation must be provided.",
                exception=BulkModulusError,
            )
        if volume_equilibrium is None:
            self.logger.error(
                "The equilibrium volume must be provided.",
                exception=BulkModulusError,
            )
        if mode is None:
            self.logger.warning(
                "The mode is not provided. The default mode 'simple' will be used.",
                exception=BulkModulusWarning,
            )
        if len(pressures_perturbation) != len(volumes_perturbation):
            self.logger.error(
                "The number of perturbation pressure and volume must be the same.",
                exception=BulkModulusError,
            )
        if pressures_perturbation.ndim == 1:
            self.logger.warning(
                "The pressure perturbation is a 1D array. The standard deviation will be set to 0.",
                exception=BulkModulusWarning,
            )
            self._pressures_perturbation_avg = pressures_perturbation
            self._pressures_perturbation_std = np.zeros_like(
                pressures_perturbation)
            pressures_perturbation = np.column_stack(
                (pressures_perturbation, np.zeros_like(pressures_perturbation)))

        elif pressures_perturbation.shape[1] == 2:
            self.logger.warning(
                "The pressure perturbation is a 2D array. The first column will be used as the average and the second column as the standard deviation.",
                exception=BulkModulusWarning,
            )
            self._pressures_perturbation_avg = pressures_perturbation[:, 0]
            self._pressures_perturbation_std = pressures_perturbation[:, 1]
        else:
            self.logger.error(
                "The pressure perturbation must be a 1D or 2D array.",
                exception=BulkModulusError,
            )

        if len(pressures_perturbation) < 2:
            self.logger.error(
                "At least two perturbation pressures are required.",
                exception=BulkModulusError,
            )
        if len(volumes_perturbation) < 2:
            self.logger.error(
                "At least two perturbation volumes are required.",
                exception=BulkModulusError,
            )

        if not np.all(np.diff(volumes_perturbation) > 0):
            self.logger.error(
                "The volume perturbation must be sorted from smallest to largest.",
                exception=BulkModulusError,
            )

        self._volumes_perturbation = volumes_perturbation
        self._volume_equilibrium = volume_equilibrium

        if mode is None:
            self._mode = self._mode_default
        elif mode == "BMEOS" or mode == "MEOS" or mode == "simple":
            self._mode = mode
        else:
            self.logger.error(
                "The mode must be either 'simple', 'MEOS' or 'BMEOS'",
                exception=BulkModulusError,
            )

        # dummy implementation
        self._bulk_modulus = 0
        self._bulk_modulus_prime = 0
        self._bulk_modulus_err = 0
        self._bulk_modulus_prime_err = 0

    def _initialize_run_simple(self):
        """
        Initializes the bulk modulus analysis.
        """
        if len(self._volumes_perturbation) != 2:
            self.logger.error(
                "The volume perturbation must have a length of 2.",
                exception=BulkModulusError,
            )
        if len(self._pressures_perturbation_avg) != 2:
            self.logger.error(
                "The pressure perturbation must have a length of 2.",
                exception=BulkModulusError,
            )

    def _two_point_stencel(self):
        """
        Calculates the finite difference using a two-point stencil.

        .. math::
            f'(x) = \frac{f(x+h) - f(x)}{h}
        Returns
        -------
        float
            the finite difference data
        """
        return (self._pressures_perturbation_avg[1] - self._pressures_perturbation_avg[0]) / (
            self._volumes_perturbation[1] - self._volumes_perturbation[0])

    def _fit_bmeos(self):
        """
        Fits the Birch-Murnaghan equation of state to the data.
        """

        initial_guesses = [200, 100, self._volume_equilibrium]
        boundary = ([0, 0, self._volume_equilibrium-0.00000001], [
                    np.inf, np.inf, self._volume_equilibrium+0.00000001])
        popt, pcov = curve_fit(
            f=_third_order_bmeos,
            xdata=self._volumes_perturbation,
            ydata=self._pressures_perturbation_avg,
            p0=initial_guesses,
            sigma=self._pressures_perturbation_std,
            absolute_sigma=True,
            bounds=boundary
        )

        b0_opt, b0_prime_opt, v0_opt = popt
        b0_err, b0_prime_err, v0_err = np.sqrt(np.diag(pcov))

        if ((self._volume_equilibrium - v0_opt) > 0.0001) or v0_err > 0.1:
            self.logger.error("The equilibrium volume {} is not close to the fitted equilibrium volume {}+-/{}.".format(
                self._volume_equilibrium, v0_opt, v0_err),
                exception=BulkModulusError,
            )

        return b0_opt, b0_prime_opt, b0_err, b0_prime_err

    def _run_bmeos(self):
        """
        Calculates the bulk modulus using the Birch-Murnaghan equation of state.

        .. math::
            B = -V \\left(\\frac{dP}{dV}\\right)_T
        """
        self._initialize_run_simple()

        self._bulk_modulus, self._bulk_modulus_prime, self._bulk_modulus_err, self._bulk_modulus_prime_err = self._fit_bmeos()

    def _calculate_bulk_modulus_simple(self):
        """
        Calculates the bulk modulus using the finite difference method.

        .. math::
            B = -V \\left(\\frac{dP}{dV}\\right)_T
        """

        finite_difference = self._two_point_stencel()
        self._bulk_modulus = -(self._volume_equilibrium * finite_difference)
    
    def _run_simple(self):
        """
        Runs the simple finite difference method.

        Returns
        -------
        List[1DNumberArray, 1DNumberArray, 1DNumberArray]
            the volume perturbation,
            the pressure perturbation average,
            the pressure perturbation standard deviation,
        """
        self._initialize_run_simple()
        self._calculate_bulk_modulus_simple()

        return [
            self._volumes_perturbation,
            self._pressures_perturbation_avg,
            self._pressures_perturbation_std
        ]

    @ timeit_in_class
    def run(self):
        """
        Runs the bulk modulus analysis
        depending on the mode:
        simple, MEOS, BMEOS.
        """
        if self._mode == "simple":
            return self._run_simple()

        if self._mode == "MEOS":
            self.logger.error(
                "Not implemented yet",
                exception=BulkModulusError,
            )

        if self._mode == "BMEOS":
            return self._run_bmeos()

        else:
            self.logger.error(
                "The mode must be either 'simple', 'MEOS' or 'BMEOS'",
                exception=BulkModulusError,
            )


        # ###########################################################################

        # # Getters

    @ property
    def bulk_modulus(self) -> float:
        """
        float: The bulk modulus.
        """
        return self._bulk_modulus

    @ property
    def bulk_modulus_prime(self) -> float:
        """
        float: The first derivative of the bulk modulus.
        """
        return self._bulk_modulus_prime

    @ property
    def bulk_modulus_err(self) -> float:
        """
        float: The error of the bulk modulus.
        """
        return self._bulk_modulus_err

    @ property
    def bulk_modulus_prime_err(self) -> float:
        """
        float: The error of the first derivative of the bulk modulus.
        """
        return self._bulk_modulus_prime_err

    @ property
    def pressures_perturbation_avg(self) -> Np1DNumberArray:
        """
        Np1DNumberArray: The average pressure perturbation.
        """
        return self._pressures_perturbation_avg

    @ property
    def pressures_perturbation_std(self) -> Np1DNumberArray:
        """
        Np1DNumberArray: The standard deviation of the pressure perturbation.
        """
        return self._pressures_perturbation_std

    @ property
    def volumes_perturbation(self) -> Np1DNumberArray:
        """
        Np1DNumberArray: The volume perturbation.
        """
        return self._volumes_perturbation

    @ property
    def volume_equilibrium(self) -> float:
        """
        float: The equilibrium volume.
        """
        return self._volume_equilibrium

    @ property
    def mode(self) -> str:
        """
        str: The mode of the analysis.
        """
        return self._mode
