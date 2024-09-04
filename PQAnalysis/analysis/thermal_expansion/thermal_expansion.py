"""
A module containing the classes for thermal expansion analysis.
"""

import logging

from beartype.typing import List

# 3rd party imports
import numpy as np

# 3rd party imports

# local absolute imports
# from PQAnalysis.config import with_progress_bar
from PQAnalysis.types import Np1DNumberArray, Np2DNumberArray
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__

from PQAnalysis.type_checking import runtime_type_checking
# local relative imports
from .exceptions import ThermalExpansionError, ThermalExpansionWarning



class ThermalExpansion:

    """
    A class to perform linear or volumetric thermal expansion coefficient analysis.

    The thermal expansion coefficient is calculated using finite differentiation.
    At the moment only a five-point stencil is supported.
    Therefore, the number of temperature points and boxeses must be 5.
    So at 5 different temperatures, the boxes data must be provided.
    The cell data should be a list of Cells objects.

    .. math::

        \\alpha = \\frac{1}{L_{2}}\\left (\\frac{\\partial L}{\\partial T}\\right )_{P}

    .. math::

        \\left (\\frac{\\partial L}{\\partial T} \\right )_{P} =  
        \\frac{\\left<L_{0}\\right> - 8 \\left <L_{1} \\right> + 
        8 \\left<L_{3}\\right> - \\left <L_{4} \\right >}{12 \\Delta T}
    -----
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        temperature_points: Np1DNumberArray | None = None,
        boxes_avg: List[Np1DNumberArray] | None = None,
        boxes_std: List[Np1DNumberArray] | None = None
    ):
        """
        Parameters
        ----------
        temperature_points : Np1DNumberArray
            the temperature points, by default None
        boxes_avg : List[Np1DNumberArray]
            the average boxes data, by default None
        boxes_std : List[Np1DNumberArray]
            the standard deviation of the boxes data, by default None
        """
        # dummy implementation
        self._temperature_step_size = 0
        self._thermal_expansions = np.zeros(4)
        self._middle_points = np.zeros(4)

        # check temperature_points have always the same step size
        if temperature_points is not None:
            self._temperature_step_size = (
                temperature_points[1] - temperature_points[0]
            )
            self._temperature_points = temperature_points
        else:
            self.logger.error(
                "Temperature points must be provided",
                exception=ThermalExpansionError
            )
        if not np.allclose(
            np.diff(temperature_points),
            self._temperature_step_size,
            rtol=1e-4
        ):
            self.logger.error(
                "Temperature points must have the same step size",
                exception=ThermalExpansionError
            )
        if len(self._temperature_points) != 5:
            self.logger.error(
                (
                    f"The number of temperature points must be 5. "
                    f"You have provided {len(self._temperature_points)} points. "
                    f"Only 5-point stencil is supported at the moment!"
                ),
                exception=ThermalExpansionError
            )

        if boxes_avg is not None:
            boxes_avg = np.array(boxes_avg)
            if np.shape(boxes_avg)[0] == 4:
                self._boxes_avg = boxes_avg
            elif np.shape(boxes_avg)[1] == 4:
                self._boxes_avg = boxes_avg.T
            else:
                self.logger.error(
                    "The boxes data must have 4 columns",
                    exception=ThermalExpansionError
                )
            if np.shape(self._boxes_avg)[1] != len(self._temperature_points):
                self.logger.error(
                    "The boxes data must have the same length as the temperature points",
                    exception=ThermalExpansionError
                )
        else:
            self.logger.error(
                "Cell data must be provided", exception=ThermalExpansionError
            )
        if boxes_std is not None:
            boxes_std = np.array(boxes_std)
            if np.shape(boxes_std)[0] == 4:
                self._boxes_std = boxes_std
            elif np.shape(boxes_std)[1] == 4:
                self._boxes_std = boxes_std.T
            else:
                self.logger.error(
                    "The boxes std data must have 4 columns",
                    exception=ThermalExpansionError
                )
            if np.shape(self._boxes_std)[1] != len(self._temperature_points):
                self.logger.error(
                    "The boxes std data must have the same length as the temperature points",
                    exception=ThermalExpansionError
                )
        else:
            self._boxes_std = np.zeros_like(self._boxes_avg)
            self.logger.warning(
                "The standard deviation of the boxes is set to zero because no data is provided",
                exception=ThermalExpansionWarning
            )

    def _initialize_run(self):
        """
        Initializes the thermal expansion analysis.
        Calculates the average boxes data and
        the standard deviation of the boxes data.
        For the thermal expansion analysis,
        the middle points of the boxes data are calculated.

        The data is stored in the following format:
        [a_avg, b_avg, c_avg, volume_avg]
        [a_std, b_std, c_std, volume_std]
        in the self._boxes_avg and self._boxes_std attributes.
        """

        middle_point = []
        middle_point.append(
            self._boxes_avg[:, np.shape(self._boxes_avg)[0] // 2]
        )
        self._middle_points = np.array(middle_point)

    def _five_point_stencel(self):
        """
        Calculates the finite difference using a five-point stencil.

        .. math::
            f'(x) = \\frac{f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)}{12h}`

        Returns
        -------
        Np1DNumberArray
            the finite difference data
        """

        finite_difference_data = (
            self._boxes_avg[:, 0] - 8 * self._boxes_avg[:, 1] +
            8 * self._boxes_avg[:, 3] - self._boxes_avg[:, 4]
        ) / (12 * self._temperature_step_size)

        return finite_difference_data

    def _calculate_thermal_expansion(self):
        """
        Calculates the thermal expansion coefficients.

        .. math::
            \\alpha = \\left (\\frac{\\partial L}{L \\partial T} \\right )_{P}
        """

        thermal_deviations = self._five_point_stencel()
        self._thermal_expansions = thermal_deviations / self._middle_points

    @timeit_in_class
    def run(self):
        """
        Runs the thermal expansion analysis.

        Returns
        -------
        List[Np1DNumberArray]
            the average boxes data, 
            the standard deviation of the boxes data and
            the thermal expansion coefficients
        """
        self._initialize_run()
        self._calculate_thermal_expansion()
        return [self._boxes_avg, self._boxes_std, self._thermal_expansions]

    ###########################################################################

    # Getters

    @property
    def temperature_points(self):
        """
        Returns
        -------
        Np1DNumberArray
            the temperature points
        """
        return self._temperature_points

    @property
    def temperature_step_size(self):
        """
        Returns
        -------
        float
            the temperature step size
        """
        return self._temperature_step_size

    @property
    def boxes_avg(self):
        """
        Returns
        -------
        Np2DNumberArray
            the average boxes data
        """
        return Np2DNumberArray(self._boxes_avg)

    @property
    def boxes_std(self):
        """
        Returns
        -------
        Np2DNumberArray
            the standard deviation of the boxes data
        """
        return Np2DNumberArray(self._boxes_std)

    @property
    def thermal_expansions(self):
        """
        Returns
        -------
        Np1DNumberArray
            the thermal expansion coefficients
        """
        return Np1DNumberArray(self._thermal_expansions[0])

    @property
    def middle_points(self):
        """
        Returns
        -------
        Np2DNumberArray
            the middle points of the boxes data
        """
        return Np2DNumberArray(self._middle_points)
