"""
A module containing the Box class.
"""

import logging

from collections import defaultdict

import numpy as np

from beartype.typing import Dict

from PQAnalysis.types import Np2DNumberArray, Np1DNumberArray
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

from .exceptions import EnergyError


class Box():

    """
    A class to store the lattice parameter data a, b, c, alpha, beta, gamma.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        a: Np1DNumberArray,
        b: Np1DNumberArray,
        c: Np1DNumberArray,
        alpha: Np1DNumberArray,
        beta: Np1DNumberArray,
        gamma: Np1DNumberArray,
        unit: str | None = None

    ) -> None:
        """
        Parameters
        ----------
        a : Np1DNumberArray
            The lattice parameter a.
        b : Np1DNumberArray
            The lattice parameter b.
        c : Np1DNumberArray
            The lattice parameter c.
        alpha : Np1DNumberArray
            The lattice parameter alpha.
        beta : Np1DNumberArray
            The lattice parameter beta.
        gamma : Np1DNumberArray
            The lattice parameter gamma.
        unit : str, optional
            The unit of the lattice parameters, by default None
        """

        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.unit = unit

    def volume(self) -> Np1DNumberArray:
        """
        Calculate the volume of the unit cell.

        Returns
        -------
        Np1DNumberArray
            The volume of the unit cell.
        """

        volume = self.a * self.b * self.c * np.sqrt(
            1 - np.cos(self.alpha)**2 - np.cos(self.beta)**2 - np.cos(self.gamma)**2 + 2 * np.cos(self.alpha) * np.cos(self.beta) * np.cos(self.gamma))

        return volume
