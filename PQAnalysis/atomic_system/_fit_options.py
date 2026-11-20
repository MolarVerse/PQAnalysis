"""
A module containing the options for fitting a system into an AtomicSystem.
"""

from dataclasses import dataclass

from PQAnalysis.types import Np1DNumberArray, PositiveInt, PositiveReal



@dataclass(frozen=True)
class FitOptions:

    """
    Options for :py:meth:`~PQAnalysis.atomic_system.AtomicSystem.iter_fit_atomic_system`.

    Attributes
    ----------
    max_iterations : PositiveInt
        The maximum number of iterations per addition, by default 100
    distance_cutoff : PositiveReal
        The distance cutoff for the fitting, by default 1.0
    max_displacement : PositiveReal | Np1DNumberArray
        The maximum displacement percentage for the fitting, by default 0.1
    rotation_angle_step : PositiveInt
        The angle step for the rotation of the system, by default 10
    """

    max_iterations: PositiveInt = 100
    distance_cutoff: PositiveReal = 1.0
    max_displacement: PositiveReal | Np1DNumberArray = 0.1
    rotation_angle_step: PositiveInt = 10
