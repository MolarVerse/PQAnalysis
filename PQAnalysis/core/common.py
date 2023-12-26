import numpy as np

from .cell import cell
from ..types import Np1DNumberArray, Np2DNumberArray, PositiveReal


def distance(pos1: Np1DNumberArray, pos2: Np1DNumberArray | Np2DNumberArray, cell: cell.Cell = cell.Cell()) -> PositiveReal | Np1DNumberArray:
    """
    Returns the distance between two positions including periodic boundary conditions.

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
