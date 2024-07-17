"""
A module containing the classes for finite differentiation analysis.
"""

from beartype.typing import List
import logging

# 3rd party imports
import numpy as np

# 3rd party imports

from tqdm.auto import tqdm

# local absolute imports
from PQAnalysis.config import with_progress_bar
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.physical_data.box import Box
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.analysis.finite_difference import FiniteDifference
# local relative imports
from .exceptions import ThermalExpansionError


class ThermalExpansion:
    """
    A class to perform linear or volumetric thermal expansion coefficient analysis.

    """

    @runtime_type_checking
    def __init__(
        self,
        temperature_points: Np1DNumberArray | None = None,
        box: List[Box] | None = None,
    ):
        """
        Parameters
        ----------
        temperature_points : Np1DNumberArray, optional
            the temperature points, by default None
        box : Box
            the box data: a, b, c, alpha, beta, gamma
        """

        self.box = box
        self.temperature_points = temperature_points
        self.volume = self.box.volume()

    def _initialize_run(self):
        """
        Initializes the thermal expansion analysis.
        """

        self.a_av = [np.average(box.a) for box in self.box]
        self.b_av = [np.average(box.b) for box in self.box]
        self.c_av = [np.average(box.c) for box in self.box]
        self.a_std = [np.std(box.a) for box in self.box]
        self.b_std = [np.std(box.b) for box in self.box]
        self.c_std = [np.std(box.c) for box in self.box]
        self.volume_av = [box.volume() for box in self.box]
        # Assuming volume() returns an array-like structure
        self.volume_std = [np.std([box.volume()]) for box in self.box]

        self.a_av = np.array(self.a_av)
        self.b_av = np.array(self.b_av)
        self.c_av = np.array(self.c_av)
        self.a_std = np.array(self.a_std)
        self.b_std = np.array(self.b_std)
        self.c_std = np.array(self.c_std)
        self.volume_av = np.array(self.volume_av)
        self.volume_std = np.array(self.volume_std)

        self._box_av = np.array(
            [self.a_av, self.b_av, self.c_av, self.volume_av])
        self._box_std = np.array(
            [self.a_std, self.b_std, self.c_std, self.volume_std])

    def _calculate_thermal_expansion(self):
        """
        Calculates the thermal expansion coefficients.
        """

        thermal_deviation_a, thermal_deviation_a_std = FiniteDifference(
            self.temperature_points, self.a_av, self.a_std).run()
        thermal_deviation_b, thermal_deviation_b_std = FiniteDifference(
            self.temperature_points, self.b_av, self.b_std).run()
        thermal_deviation_c, thermal_deviation_c_std = FiniteDifference(
            self.temperature_points, self.c_av, self.c_std).run()
        thermal_deviation_volume, thermal_deviation_volume_std = FiniteDifference(
            self.temperature_points, self.volume_av, self.volume_std).run()

        a_middle_point = self.a_av[len(self.a_av) // 2]
        b_middle_point = self.b_av[len(self.b_av) // 2]
        c_middle_point = self.c_av[len(self.c_av) // 2]
        volume_middle_point = self.volume_av[len(self.volume_av) // 2]

        thermal_expansion_a = thermal_deviation_a / a_middle_point
        thermal_expansion_b = thermal_deviation_b / b_middle_point
        thermal_expansion_c = thermal_deviation_c / c_middle_point
        thermal_expansion_volume = thermal_deviation_volume / volume_middle_point

        self._thermal_expansion = np.array(
            [thermal_expansion_a, thermal_expansion_b, thermal_expansion_c, thermal_expansion_volume])

    def run(self):
        """
    Runs the thermal expansion analysis.

    Returns
    -------
    Np1DNumberArray
        the linear or volumetric thermal expansion coefficients
    """
        self._initialize_run()
        self._calculate_thermal_expansion()
        return [self._box_av, self._box_std, self._thermal_expansion]
