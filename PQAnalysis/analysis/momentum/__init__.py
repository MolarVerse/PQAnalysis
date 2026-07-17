"""
A package containing classes and functions to check the total linear
momentum of velocity trajectories.

Classes
-------
:py:class:`~PQAnalysis.analysis.momentum.momentum.Momentum`
    A class to calculate the norm of the total linear momentum of a
    selection of atoms for every frame of a velocity trajectory.
:py:class:`~PQAnalysis.analysis.momentum.momentum_output_file_writer.MomentumDataWriter`
    A class to write momentum data to output files.

Functions
---------
:py:func:`~PQAnalysis.analysis.momentum.api.check_momentum`
    A function to calculate and write the total linear momentum norm
    per frame of a velocity trajectory.
"""

from .api import check_momentum
from .momentum import Momentum
from .momentum_output_file_writer import MomentumDataWriter

__all__ = [
    "Momentum",
    "MomentumDataWriter",
    "check_momentum",
]
