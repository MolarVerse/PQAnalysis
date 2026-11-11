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

#: Elements from this atomic number on (X, Q, Cav, Sup) are placeholders
#: with a made-up mass, not real atoms.
PLACEHOLDER_ATOMIC_NUMBER = 999



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



def group_mass_weighted_mean(
    values: np.ndarray,
    masses: np.ndarray,
    group_ids: np.ndarray,
) -> np.ndarray:
    """
    Calculates the mass weighted mean of vectors for every group of atoms.

    This is the velocity of the center of mass for velocities. Unlike the
    positions, vectors need no periodic images.

    Parameters
    ----------
    values : np.ndarray
        The vector of every atom with shape (n_atoms, 3).
    masses : np.ndarray
        The masses of the atoms with shape (n_atoms,).
    group_ids : np.ndarray
        The group of every atom with shape (n_atoms,).

    Returns
    -------
    np.ndarray
        The mass weighted means with shape (n_groups, 3), one row per
        distinct group id in ascending order of the ids.
    """

    values = np.asarray(values, dtype=float)
    masses = np.asarray(masses, dtype=float)
    _, inverse = np.unique(np.asarray(group_ids), return_inverse=True)
    inverse = inverse.reshape(-1)
    total_mass = np.bincount(inverse, weights=masses)

    return np.stack(
        [
            np.bincount(inverse, weights=values[:, k] * masses)
            for k in range(3)
        ],
        axis=1,
    ) / total_mass[:, None]



def group_sums(values: np.ndarray, group_ids: np.ndarray) -> np.ndarray:
    """
    Calculates the sum of the values of every group of atoms.

    This is the force on the center of mass for forces and the total charge
    of a group for charges.

    Parameters
    ----------
    values : np.ndarray
        The value of every atom with shape (n_atoms,) or (n_atoms, 3).
    group_ids : np.ndarray
        The group of every atom with shape (n_atoms,).

    Returns
    -------
    np.ndarray
        The sums with shape (n_groups,) or (n_groups, 3), one row per
        distinct group id in ascending order of the ids.
    """

    values = np.asarray(values, dtype=float)
    _, inverse = np.unique(np.asarray(group_ids), return_inverse=True)
    inverse = inverse.reshape(-1)

    if values.ndim == 1:
        return np.bincount(inverse, weights=values)

    return np.stack(
        [np.bincount(inverse, weights=values[:, k]) for k in range(3)],
        axis=1,
    )



# The arrays are the natural arguments of the function.
# pylint: disable-next=too-many-arguments,too-many-positional-arguments
def center_of_mass_properties(
    masses: np.ndarray,
    group_ids: np.ndarray,
    cell: Cell,
    positions: np.ndarray | None = None,
    velocities: np.ndarray | None = None,
    forces: np.ndarray | None = None,
    charges: np.ndarray | None = None,
    check_extent: bool = False,
) -> dict:
    """
    Calculates the properties of the center of mass of every group of atoms
    for all the data that is given.

    The center of mass has the mass weighted mean of the velocities, the
    sum of the forces on the atoms as force, and the sum of the charges as
    charge. Data that is None or empty is left out.

    Parameters
    ----------
    masses : np.ndarray
        The masses of the atoms with shape (n_atoms,).
    group_ids : np.ndarray
        The group of every atom with shape (n_atoms,).
    cell : Cell
        The cell of the system used for the periodic boundary conditions.
    positions : np.ndarray, optional
        The positions with shape (n_atoms, 3).
    velocities : np.ndarray, optional
        The velocities with shape (n_atoms, 3).
    forces : np.ndarray, optional
        The forces with shape (n_atoms, 3).
    charges : np.ndarray, optional
        The charges with shape (n_atoms,).
    check_extent : bool, optional
        See :py:func:`center_of_mass_groups`, by default False.

    Returns
    -------
    dict
        The arrays of the data that was given, with the keys ``pos``,
        ``vel``, ``forces`` and ``charges`` and one row per distinct group
        id in ascending order of the ids.

    Raises
    ------
    AtomicSystemPositionsError
        If given data does not have the same number of atoms as the masses.
    AtomicSystemError
        For the errors of :py:func:`center_of_mass_groups`.
    """

    given = {
        "pos": positions,
        "vel": velocities,
        "forces": forces,
        "charges": charges,
    }
    given = {
        key: np.asarray(value, dtype=float)
        for key, value in given.items()
        if value is not None and np.size(value) > 0
    }

    names = {
        "pos": "positions",
        "vel": "velocities",
        "forces": "forces",
        "charges": "charges",
    }

    for key, value in given.items():
        if len(value) != len(masses):
            logger.error(
                (
                    f"The {names[key]} of {len(value)} atoms do not fit to "
                    f"the {len(masses)} masses."
                ),
                exception=AtomicSystemPositionsError,
            )

    if np.any(group_sums(masses, group_ids) <= 0.0):
        logger.error(
            "The total mass of every group must be positive.",
            exception=AtomicSystemError,
        )

    result = {}

    if "pos" in given:
        result["pos"] = center_of_mass_groups(
            given["pos"], masses, group_ids, cell, check_extent
        )
    if "vel" in given:
        result["vel"] = group_mass_weighted_mean(
            given["vel"], masses, group_ids
        )
    if "forces" in given:
        result["forces"] = group_sums(given["forces"], group_ids)
    if "charges" in given:
        result["charges"] = group_sums(given["charges"], group_ids)

    return result



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



def placeholder_atom_names(atoms, indices) -> list:
    """
    Returns the names of the atoms with a placeholder mass.

    Atoms of the elements X, Q, Cav and Sup have a made-up mass of 999 or
    more. They dominate the center of mass of their molecule.

    Parameters
    ----------
    atoms : list of Atom
        The atoms of the system.
    indices : np.ndarray
        The indices of the atoms to check.

    Returns
    -------
    list of str
        The sorted, distinct names of the atoms with a placeholder mass.
    """

    return sorted({
        atoms[index].name
        for index in indices
        if (atoms[index].atomic_number or 0) >= PLACEHOLDER_ATOMIC_NUMBER
    })
