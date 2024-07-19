"""
A module containing the classes for finite differentiation analysis.
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
    Therefore, the number of temperature points and boxes must be 5.
    So at 5 different temperatures, the box data must be provided.
    The box data should be a list of Box objects.

    .. math::
        \alpha = \frac{1}{L_{2}}frac{\Delta L}{\Delta T}

    .. math::
        (\frac{\Delta L}{\Delta T})_P = 
            \frac{<L_{0}> - 8 L_{1}> + 8 <L_{3}> - <L_{4}>}{12 \Delta T}
    -----
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        temperature_points: Np1DNumberArray | None = None,
        box: List[Box] | None = None,
    ):
        """
        Parameters
        ----------
        temperature_points : Np1DNumberArray
            the temperature points, by default None
        box : List[Box]
            the box data: a, b, c, alpha, beta, gamma
        """
        # dummy implementation
        self._temperature_step_size = 0
        self._box_avg = np.zeros(4)
        self._box_std = np.zeros(4)
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
        if np.all(np.diff(self._temperature_points) != self._temperature_step_size):
            self.logger.error(
                "Temperature points must have the same step size",
                exception=ThermalExpansionError
            )

        if box is not None:
            self._box = box
        else:
            self.logger.error(
                "Box data must be provided",
                exception=ThermalExpansionError
            )
        if len(self._temperature_points) != 5:
            self.logger.error(
                (
                    "The number of temperature points must be 5"
                    f"You have provided {len(self._temperature_points)} points"
                    "Only 5-point stencil is supported at the moment!"
                ),
                exception=ThermalExpansionError
            )
        if len(self._temperature_points) != len(self._box):
            self.logger.error(
                "Temperature points and box data must have the same length",
                exception=ThermalExpansionError
            )

    def _initialize_run(self):
        """
        Initializes the thermal expansion analysis.
        Calculates the average box data and 
        the standard deviation of the box data.
        For the thermal expansion analysis, 
        the middle points of the box data are calculated.

        The data is stored in the following format:
        [a_avg, b_avg, c_avg, volume_avg]
        [a_std, b_std, c_std, volume_std]
        in the self._box_avg and self._box_std attributes.
        """

        a_avg = np.array([np.average(box.a) for box in self._box])
        b_avg = np.array([np.average(box.b) for box in self._box])
        c_avg = np.array([np.average(box.c) for box in self._box])
        a_std = np.array([np.std(box.a) for box in self._box])
        b_std = np.array([np.std(box.b) for box in self._box])
        c_std = np.array([np.std(box.c) for box in self._box])
        volume_avg = np.array([np.average(box.volume()) for box in self._box])
        volume_std = np.array([np.std(box.volume()) for box in self._box])
        middle_point = []
        middle_point.append(a_avg[len(a_avg) // 2])
        middle_point.append(b_avg[len(b_avg) // 2])
        middle_point.append(c_avg[len(c_avg) // 2])
        middle_point.append(volume_avg[len(volume_avg) // 2])

        self._middle_points = np.array(middle_point)

        self._box_avg = np.array(
            [a_avg, b_avg, c_avg, volume_avg])
        self._box_std = np.array(
            [a_std, b_std, c_std, volume_std])

    def _five_point_stencel(self):
        """
        Calculates the finite difference using a five-point stencil.

        .. math::
            f'(x) = \frac{f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)}{12h}`

        Returns
        -------
        Np1DNumberArray
            the finite difference data
        """

        finite_difference_data = (self._box_avg[:, 0]-8*self._box_avg[:, 1]+8 *
                                  self._box_avg[:, 3]-self._box_avg[:, 4]) / (12 * self._temperature_step_size)

        return finite_difference_data

    def _calculate_thermal_expansion(self):
        """
        Calculates the thermal expansion coefficients.
        :math:`\\alpha = \\frac{\\Delta L}{L \\Delta T}`
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
            the average box data, 
            the standard deviation of the box data and
            the thermal expansion coefficients
        """
        self._initialize_run()
        self._calculate_thermal_expansion()
        return [self._box_avg, self._box_std, self._thermal_expansions]

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
    def box(self):
        """
        Returns
        -------
        List[Box]
            the box data
        """
        return self._box

    @property
    def box_avg(self):
        """
        Returns
        -------
        Np1DNumberArray
            the average box data
        """
        return self._box_avg

    @property
    def box_std(self):
        """
        Returns
        -------
        Np1DNumberArray
            the standard deviation of the box data
        """
        return self._box_std

    @property
    def thermal_expansions(self):
        """
        Returns
        -------
        Np1DNumberArray
            the thermal expansion coefficients
        """
        return self._thermal_expansions

    @property
    def middle_points(self):
        """
        Returns
        -------
        Np1DNumberArray
            the middle points of the box data
        """
        return self._middle_points
