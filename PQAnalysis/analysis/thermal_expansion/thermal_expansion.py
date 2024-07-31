"""
A module containing the classes for thermal expansion analysis.
"""

from beartype.typing import List
import logging

# 3rd party imports
import numpy as np

# 3rd party imports

# from tqdm.auto import tqdm

# local absolute imports
# from PQAnalysis.config import with_progress_bar
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.physical_data.box import Box

from PQAnalysis.type_checking import runtime_type_checking
# local relative imports
from .exceptions import ThermalExpansionError


class ThermalExpansion:
    """
    A class to perform linear or volumetric thermal expansion coefficient analysis.

    The thermal expansion coefficient is calculated using finite differentiation.
    At the moment only a five-point stencil is supported.
    Therefore, the number of temperature points and boxeses must be 5.
    So at 5 different temperatures, the boxes data must be provided.
    The boxes data should be a list of Box objects.

    .. math::

        \\alpha = \\frac{1}{L_{2}}\\left (\\frac{\\partial L}{\\partial T}\\right )_{P}

    .. math::

        \\left (\\frac{\\partial L}{\\partial T} \\right )_{P} =  \\frac{\\left<L_{0}\\right> - 8 \\left <L_{1} \\right> + 8 \\left<L_{3}\\right> - \\left <L_{4} \\right >}{12 \\Delta T}
    -----
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        temperature_points: Np1DNumberArray | None = None,
        boxes: List[Box] | None = None,
    ):
        """
        Parameters
        ----------
        temperature_points : Np1DNumberArray
            the temperature points, by default None
        boxes : List[Box]
            the boxes data: a, b, c, alpha, beta, gamma
        """
        # dummy implementation
        self._temperature_step_size = 0
        self._boxes_avg = np.zeros(4)
        self._boxes_std = np.zeros(4)
        self._thermal_expansions = np.zeros(4)
        self._middle_points = np.zeros(4)

        # check temperature_points have always the same step size
        if temperature_points is not None:
            self._temperature_step_size = temperature_points[1] - \
                temperature_points[0]
            self._temperature_points = temperature_points
        else:
            self.logger.error(
                "Temperature points must be provided",
                exception=ThermalExpansionError
            )
        if not np.allclose(np.diff(temperature_points), self._temperature_step_size, rtol=1e-4):
            self.logger.error(
                "Temperature points must have the same step size",
                exception=ThermalExpansionError
            )

        if boxes is not None:
            self._boxes = boxes
        else:
            self.logger.error(
                "Box data must be provided",
                exception=ThermalExpansionError
            )
        if len(self._temperature_points) != 5:
            self.logger.error(
                (
                    "The number of temperature points must be 5. "
                    f"You have provided {
                        len(self._temperature_points)} points. "
                    "Only 5-point stencil is supported at the moment!"
                ),
                exception=ThermalExpansionError
            )
        if len(self._temperature_points) != len(self._boxes):
            self.logger.error(
                "Temperature points and boxes data must have the same length",
                exception=ThermalExpansionError
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

        a_avg = np.array([np.average(boxes.a) for boxes in self._boxes])
        b_avg = np.array([np.average(boxes.b) for boxes in self._boxes])
        c_avg = np.array([np.average(boxes.c) for boxes in self._boxes])
        a_std = np.array([np.std(boxes.a) for boxes in self._boxes])
        b_std = np.array([np.std(boxes.b) for boxes in self._boxes])
        c_std = np.array([np.std(boxes.c) for boxes in self._boxes])
        volume_avg = np.array([np.average(boxes.volume())
                               for boxes in self._boxes])
        volume_std = np.array([np.std(boxes.volume())
                               for boxes in self._boxes])
        middle_point = []
        middle_point.append(a_avg[len(a_avg) // 2])
        middle_point.append(b_avg[len(b_avg) // 2])
        middle_point.append(c_avg[len(c_avg) // 2])
        middle_point.append(volume_avg[len(volume_avg) // 2])

        self._middle_points = np.array(middle_point)

        self._boxes_avg = np.array(
            [a_avg, b_avg, c_avg, volume_avg])
        self._boxes_std = np.array(
            [a_std, b_std, c_std, volume_std])

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

        finite_difference_data = (self._boxes_avg[:, 0]-8*self._boxes_avg[:, 1]+8 *
                                  self._boxes_avg[:, 3]-self._boxes_avg[:, 4]) / (12 * self._temperature_step_size)

        return finite_difference_data

    def _calculate_thermal_expansion(self):
        """
        Calculates the thermal expansion coefficients.

        .. math::
            \\alpha = \\left (\\frac{\\partial L}{L \\partial T} \\right )_{P}
        """

        thermal_deviations = self._five_point_stencel()
        self._thermal_expansions = thermal_deviations / self._middle_points

    @ timeit_in_class
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

    @ property
    def temperature_points(self):
        """
        Returns
        -------
        Np1DNumberArray
            the temperature points
        """
        return self._temperature_points

    @ property
    def temperature_step_size(self):
        """
        Returns
        -------
        float
            the temperature step size
        """
        return self._temperature_step_size

    @ property
    def boxes(self):
        """
        Returns
        -------
        List[Box]
            the boxes data
        """
        return self._boxes

    @ property
    def boxes_avg(self):
        """
        Returns
        -------
        Np1DNumberArray
            the average boxes data
        """
        return self._boxes_avg

    @ property
    def boxes_std(self):
        """
        Returns
        -------
        Np1DNumberArray
            the standard deviation of the boxes data
        """
        return self._boxes_std

    @ property
    def thermal_expansions(self):
        """
        Returns
        -------
        Np1DNumberArray
            the thermal expansion coefficients
        """
        return self._thermal_expansions

    @ property
    def middle_points(self):
        """
        Returns
        -------
        Np1DNumberArray
            the middle points of the boxes data
        """
        return self._middle_points
