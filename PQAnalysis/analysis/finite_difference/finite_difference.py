"""
A module containing the classes for finite differentiation analysis.
"""


import logging

# 3rd party imports
import numpy as np

# 3rd party imports
from beartype.typing import Tuple
from tqdm.auto import tqdm

# local absolute imports
from PQAnalysis.config import with_progress_bar
from PQAnalysis.types import Np1DNumberArray, PositiveInt, PositiveReal
from PQAnalysis.core import distance, Cells
from PQAnalysis.traj import Trajectory, check_trajectory_pbc, check_trajectory_vacuum
from PQAnalysis.topology import Selection, SelectionCompatible
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import TrajectoryReader
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

# local relative imports
from .exceptions import FiniteDifferenceError


class FiniteDifference:
    """
    A class to perform finite difference analysis.

    The finite difference Analysis calculates the derivative of a physical property by using finite difference.

    Finite difference is a numerical method to approximate the derivative of a function.
    The derivative of a function f(x) at a point x is approximated by the following formula:
    math::
        f'(x) = (f(x + h) - f(x - h)) / (2 * h)

        where h is a small number.  

    The finite difference analysis calculates the derivative of a physical property at a given point.

    """

    @runtime_type_checking
    def __init__(
        self,
        temperature_points: Np1DNumberArray | None = None,
        finite_difference_points: Tuple[Np1DNumberArray] | Np1DNumberArray | None = None,
        std_points: Tuple[Np1DNumberArray] | Np1DNumberArray | None = None,
    ):
        """
        Parameters
        ----------
        temperature_points : Np1DNumberArray, optional
            the temperature points, by default None
        finite_difference_points : Tuple[Np1DNumberArray] | Np1DNumberArray | None
            the points at which the finite difference is calculated
        std_points : Tuple[Np1DNumberArray] | Np1DNumberArray | None, optional
            the standard deviation points, by default None
        """

        if not isinstance(finite_difference_points, Tuple):
            finite_difference_points = (finite_difference_points,)
        if len(finite_difference_points) != 1 or finite_difference_points == None:
            raise FiniteDifferenceError(
                "For the finite difference analysis at least one list of points has to be provided."
            )
        if not all(
            [
                finite_difference_points[0].shape == points.shape
                for points in finite_difference_points
            ]
        ):
            raise FiniteDifferenceError(
                "All list of points must have the same shape."
            )
        # check the size of the points
        if finite_difference_points[0].shape[0] != 5:
            raise FiniteDifferenceError(
                "The finite difference analyis is implemented only for 5 points."
            )
        if not isinstance(temperature_points, np.ndarray):
            raise FiniteDifferenceError(
                "The temperature points must be a numpy array."
            )

        if std_points is not None:
            if not isinstance(std_points, Tuple):
                std_points = (std_points,)
            if len(std_points) != 1 or std_points == None:
                raise FiniteDifferenceError(
                    "For the standard deviation analysis at least one list of points has to be provided."
                )
            if not all(
                [
                    std_points[0].shape == points.shape for points in std_points
                ]
            ):
                raise FiniteDifferenceError(
                    "All list of points must have the same shape."
                )
            if std_points[0].shape[0] != 5:
                raise FiniteDifferenceError(
                    "The standard deviation analyis is implemented only for 5 points."
                )
        self._finite_difference_points = finite_difference_points
        self._temperature_points = temperature_points

    @property
    def finite_difference_points(self) -> Tuple[Np1DNumberArray]:
        """
        Returns the finite difference points.
        """
        return self._finite_difference_points

    def __str__(self) -> str:
        """
        Returns a string representation of the finite difference points.
        """
        return "\n".join(
            [
                f"Finite Difference Points: {points}"
                for points in self._finite_difference_points
            ]
        )

    @timeit_in_class
    def run(self) -> [Np1DNumberArray, Np1DNumberArray]:
        """
        Run the finite difference analysis.

        Returns
        -------
        Np1DNumberArray
            the finite difference data
        """
        finite_difference_data = self._finite_difference()
        finite_difference_std = self._finite_difference_std()

        return [finite_difference_data, finite_difference_std]

    def _finite_difference(self) -> Np1DNumberArray:
        """
        Calculate the finite difference.

        Returns
        -------
        Np1DNumberArray
            the finite difference data
        """
        finite_difference_data = np.zeros(len(self._finite_difference_points))

        for i, point in enumerate(self._finite_difference_points):
            finite_difference_data[i] = (
                point[0] - 8 * point[1] + 8 * point[3] - point[4]
            ) / 12
        return finite_difference_data

    def _finite_difference_std(self) -> Np1DNumberArray:
        """
        Calculate the standard deviation of the finite difference.

        Returns
        -------
        Np1DNumberArray
            the standard deviation of the finite difference data
        """
        if self._std_points is None:
            return None

        finite_difference_std = np.zeros(len(self._std_points))

        # TODO: Implement the standard deviation calculation
        return finite_difference_std
