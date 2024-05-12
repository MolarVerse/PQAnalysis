"""
A module containing api functions to handle trajectories and frame objects.
"""

from PQAnalysis.core import Cells
from PQAnalysis.type_checking import runtime_type_checking



@runtime_type_checking
def check_trajectory_pbc(cells: Cells) -> bool:
    """
    Checks no cell in the trajectory is Cell() i.e. 
    checks if the trajectory is never in vacuum.

    Parameters
    ----------
    cells : Cells
        The list of cells of the trajectory.

    Returns
    -------
    bool
        False if one cell of the trajectory is Cell(), True otherwise.
    """

    if len(cells) == 0:
        return False

    return all(not cell.is_vacuum for cell in cells)



@runtime_type_checking
def check_trajectory_vacuum(cells: Cells) -> bool:
    """
    Checks if all cells of the trajectory are in vacuum i.e. cell = Cell().

    Parameters
    ----------
    cells : Cells
        The list of cells of the trajectory.

    Returns
    -------
    bool
        True if all cells of the trajectory are in vacuum, False otherwise.
    """

    return not any(not cell.is_vacuum for cell in cells)
