"""
A module containing all api functions for the PQAnalysis.core package.

All of this functions are provided within the namespace of the PQAnalysis.core package.
"""

import numpy as np

from PQAnalysis.core import cell as _cell
from PQAnalysis.types import Np1DNumberArray, Np2DNumberArray, PositiveReal
from PQAnalysis.type_checking import runtime_type_checking



@runtime_type_checking
def distance(
    pos1: Np1DNumberArray | Np2DNumberArray,
    pos2: Np1DNumberArray | Np2DNumberArray,
    cell: _cell.Cell = _cell.Cell(),
    **kwargs  # pylint: disable=unused-argument
) -> PositiveReal | Np1DNumberArray | Np2DNumberArray:
    """
    Returns the distances between all combinations of two position arrays.

    Notes
    -----
    In more detail this function returns the distances between all combinations
    of two positions in the position arrays 'pos1' and 'pos2'. If pos1 and pos2
    are both 1D arrays, the distance between the two positions is returned. 
    If pos1 is a 1D array and pos2 is a 2D array, the distances between the 
    position in pos1 and all positions in pos2 are returned. If pos1 and pos2
    are both 2D arrays, the distances between all combinations of positions in
    pos1 and pos2 are returned.

    The only requirement for the position arrays is that the last dimension has
    to be 3, which represents the x, y and z coordinates of the position(s).


    Examples
    --------

    >>> from PQAnalysis.core import distance, Cell
    >>> import numpy as np

    >>> pos1 = np.array([0, 0, 0])
    >>> pos2 = np.array([1, 1, 1])
    >>> distance(pos1, pos2)
    array([[1.73205081]])

    >>> pos1 = np.array([0, 0, 0])
    >>> pos2 = np.array([[1, 1, 1], [0.5, 0.5, 0.5]])
    >>> distance(pos1, pos2)
    array([[1.73205081, 0.8660254 ]])

    >>> pos1 = np.array([0, 0, 0])
    >>> pos2 = np.array([[1, 1, 1], [0.5, 0.5, 0.5]])
    >>> distance(pos1, pos2, cell=Cell(0.7, 0.7, 0.7))
    array([[0.519615, 0.3464102]])

    Parameters
    ----------
    pos1 : Np1DNumberArray | Np2DNumberArray
        The first position.
    pos2 : Np1DNumberArray | Np2DNumberArray
        The second position.
    cell : Cell, optional
        The unit cell of the system. Default is Cell().
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    PositiveReal | Np1DNumberArray | Np2DNumberArray
        The distance(s) between the two position(s) (arrays).
    """

    pos1 = np.atleast_2d(pos1)

    delta_pos = pos2 - pos1[:, None]

    delta_pos = cell.image(delta_pos)

    return np.linalg.norm(delta_pos, axis=-1)
