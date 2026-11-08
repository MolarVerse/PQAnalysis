"""
A module containing the center of mass of groups of atoms.

The function of this module works on plain arrays, so that it can be used
without constructing an :py:class:`~PQAnalysis.atomic_system.atomic_system.AtomicSystem`
for every group (molecule, residue, block of atoms) of every frame.
"""

import logging

import numpy as np

from PQAnalysis import __package_name__
from PQAnalysis.core import Cell
from PQAnalysis.utils.custom_logging import setup_logger

from .exceptions import AtomicSystemError, AtomicSystemPositionsError

logger = logging.getLogger(__package_name__).getChild("center_of_mass_groups")
logger = setup_logger(logger)



def center_of_mass_groups(
    positions: np.ndarray,
    masses: np.ndarray,
    group_ids: np.ndarray,
    cell: Cell,
    check_extent: bool = False,
) -> np.ndarray:
    """
    Calculates the center of mass of every group of atoms.

    Every atom is imaged relative to the first atom of its group (the first
    one in the order of the arrays), so that groups which are split by the
    periodic boundary get the center of mass of the whole group. The centers
    of mass are imaged into the cell. A group therefore has to be smaller than
    half of the box in every direction, otherwise its atoms cannot be told
    apart from their periodic images and the center of mass is wrong. This
    can only be partly checked (see check_extent), and not at all by default,
    because a center of mass of a whole periodic system is also computed
    with this function.

    Parameters
    ----------
    positions : np.ndarray
        The positions of the atoms in Angstrom with shape (n_atoms, 3).
    masses : np.ndarray
        The masses of the atoms with shape (n_atoms,).
    group_ids : np.ndarray
        The group of every atom with shape (n_atoms,). The ids have to be
        integers, but they do not have to be sorted, contiguous or start at 0.
    cell : Cell
        The cell of the system used for the periodic boundary conditions.
    check_extent : bool, optional
        If True, a group whose atoms span half of the box or more (in
        fractional coordinates, along any cell vector) raises an error. The
        atoms are compared after they are imaged around the first atom of the
        group, so the check finds groups with atoms on both sides of the
        first atom, like a 6.4 Angstrom molecule in a 10 Angstrom box, but it
        cannot find a pair of atoms that are farther apart than half the
        box, since this is not distinguishable from their periodic image.
        A group smaller than half the box never raises. By default False.

    Returns
    -------
    np.ndarray
        The centers of mass in Angstrom with shape (n_groups, 3), one row
        per distinct group id in ascending order of the ids.

    Raises
    ------
    AtomicSystemPositionsError
        If positions, masses and group_ids do not describe the same number
        of atoms.
    AtomicSystemError
        If the positions do not have the shape (n_atoms, 3), if the total
        mass of a group is not positive or, with check_extent, if a group
        spans half of the box or more.
    """

    positions = np.asarray(positions, dtype=float)
    masses = np.asarray(masses, dtype=float)
    group_ids = np.asarray(group_ids)

    if positions.ndim != 2 or positions.shape[1] != 3:
        logger.error(
            "The positions must have the shape (n_atoms, 3), "
            f"but have the shape {positions.shape}.",
            exception=AtomicSystemError,
        )

    if not len(positions) == len(masses) == len(group_ids):
        logger.error(
            (
                "The number of positions, masses and group ids must be equal, "
                f"but are {len(positions)}, {len(masses)} and "
                f"{len(group_ids)}."
            ),
            exception=AtomicSystemPositionsError,
        )

    unique_ids, first_atom, inverse = np.unique(
        group_ids, return_index=True, return_inverse=True
    )
    inverse = inverse.reshape(-1)

    total_mass = np.bincount(inverse, weights=masses)

    if np.any(total_mass <= 0.0):
        logger.error(
            "The total mass of every group must be positive.",
            exception=AtomicSystemError,
        )

    reference = positions[first_atom]
    relative = cell.image(positions - reference[inverse])

    if check_extent:
        _check_extent(relative, inverse, unique_ids, cell)

    weighted = np.stack(
        [
            np.bincount(inverse, weights=relative[:, k] * masses)
            for k in range(3)
        ],
        axis=1,
    )

    return cell.image(reference + weighted / total_mass[:, None])



def _check_extent(
    relative: np.ndarray,
    inverse: np.ndarray,
    group_ids: np.ndarray,
    cell: Cell,
) -> None:
    """
    Raises an error for a group that spans half of the box or more.

    The positions are given relative to the first atom of the group and
    imaged, so the extent of a group is the range of its fractional
    coordinates along every cell vector.
    """

    if len(relative) == 0:
        return

    fractional = relative @ cell.inverse_box_matrix.T

    order = np.argsort(inverse, kind="stable")
    starts = np.searchsorted(inverse[order], np.arange(len(group_ids)))
    fractional = fractional[order]

    extent = (
        np.maximum.reduceat(fractional, starts, axis=0) -
        np.minimum.reduceat(fractional, starts, axis=0)
    )
    too_large = np.flatnonzero(np.any(extent >= 0.5, axis=1))

    if len(too_large) > 0:
        logger.error(
            (
                f"The group {group_ids[too_large[0]]} spans half of the box "
                "or more, its center of mass cannot be determined with "
                "periodic boundary conditions "
                f"({len(too_large)} of {len(group_ids)} groups are affected)."
            ),
            exception=AtomicSystemError,
        )
