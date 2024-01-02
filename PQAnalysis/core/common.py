"""
A module containing common functions for the PQAnalysis.core package.

All of this functions are provided within the namespace of the PQAnalysis.core package.

Functions
---------
distance
    Returns the distance between one position and one or more positions including periodic boundary conditions.
    
    
"""

import numpy as np

from . import cell
from ..types import Np1DNumberArray, Np2DNumberArray, PositiveReal


def distance(pos1: Np1DNumberArray, pos2: Np1DNumberArray | Np2DNumberArray, cell: cell.Cell = cell.Cell()) -> PositiveReal | Np1DNumberArray:
    """
    Returns the distance between one position and one or more positions including periodic boundary conditions.

    Examples
    --------

    >>> from PQAnalysis.core import distance
    >>> import numpy as np
    >>> pos1 = np.array([0, 0, 0])
    >>> pos2 = np.array([1, 1, 1])
    >>> distance(pos1, pos2)
    1.7320508075688772

    >>> pos1 = np.array([0, 0, 0])
    >>> pos2 = np.array([[1, 1, 1], [0.5, 0.5, 0.5]])
    >>> distance(pos1, pos2)
    array([1.73205081, 0.8660254 ])

    >>> pos1 = np.array([0, 0, 0])
    >>> pos2 = np.array([[1, 1, 1], [0.5, 0.5, 0.5]])
    >>> distance(pos1, pos2, cell=cell.Cell(0.7, 0.7, 0.7))
    array([0.519615, 0.3464102])

    Parameters
    ----------
    pos1 : Np1DNumberArray
        The first position.
    pos2 : Np1DNumberArray | Np2DNumberArray
        The second position.
    cell : Cell, optional
        The unit cell of the system. Default is Cell().

    Returns
    -------
    PositiveReal | Np1DNumberArray
        The distance(s) between the two positions.
    """

    delta_pos = pos2-pos1

    delta_pos = cell.image(delta_pos)

    return np.linalg.norm(delta_pos, axis=-1)
