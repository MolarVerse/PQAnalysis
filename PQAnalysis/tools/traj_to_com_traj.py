"""
A module containing the tool to compute a center of mass trajectory for a given selection.
"""

import logging

import numpy as np

from PQAnalysis import __package_name__
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.atomic_system._center_of_mass import center_of_mass_groups
from PQAnalysis.core import Atom
from PQAnalysis.exceptions import PQValueError
from PQAnalysis.topology import Topology
from PQAnalysis.traj import Trajectory
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.utils.custom_logging import setup_logger

logger = logging.getLogger(__package_name__).getChild("traj_to_com_traj")
logger = setup_logger(logger)



# TODO: add atom to element mapper if atom name not element names
# TODO rethink concept of selection/getitem/slices
@runtime_type_checking
def traj_to_com_traj(
    trajectory: Trajectory,
    selection=None,
    group=None,
    by_residue: bool = False,
    topology: Topology | None = None,
) -> Trajectory:
    """
    Function that computes the center of mass trajectory for a given selection.

    With the group parameter, the trajectory will be group according to the given group.
    For example if group = 2, the trajectory will be grouped in pairs of 2 and the center
    of mass of each pair will be computed.

    With by_residue, the atoms are grouped by the residues (molecules) of the
    topology, i.e. by the molecules defined in the moldescriptor file. This
    works for molecules of different sizes and types. The center of mass of
    a molecule is computed from its selected atoms; molecules without a
    selected atom are left out. The result has one atom per molecule, ordered
    by the residue number, named by the names of the selected atoms.

    Parameters
    ----------
    trajectory : Trajectory
        The trajectory to compute the center of mass trajectory for.
    selection : Atom | int | slice | Np1DIntArray, optional
        The selection to compute the center of mass trajectory for.
        If None, the selection is all atoms.
    group : int, optional
        The group to compute the center of mass trajectory for.
        If None, the group is all atoms.
    by_residue : bool, optional
        If True, the centers of mass of the residues of the topology are
        computed instead of the ones of groups of consecutive atoms.
        Cannot be combined with group. By default False.
    topology : Topology, optional
        The topology with the residues to use with by_residue. If None, the
        topology of the trajectory is used. It has to contain residues, e.g.
        by reading a restart file together with the moldescriptor file.

    Returns
    -------
    Trajectory
        The center of mass trajectory.

    Raises
    ------
    PQValueError
        If group and by_residue are combined, if by_residue is used with a
        topology without residues or if the number of atoms of the topology
        and of a frame differ.
    """

    if by_residue and group is not None:
        logger.error(
            "The parameters group and by_residue cannot be combined.",
            exception=PQValueError,
        )

    if len(trajectory) == 0:
        return Trajectory()

    if by_residue:
        return _traj_to_com_traj_by_residue(trajectory, selection, topology)

    com_traj = Trajectory()
    for frame in trajectory:
        selected = frame if selection is None else frame[selection]
        com_traj.append(selected.compute_com_atomic_system(group=group))

    return com_traj


def _selected_indices(atoms, selection) -> np.ndarray:
    """Returns the indices of the atoms, as AtomicSystem.__getitem__ does."""

    if selection is None:
        return np.arange(len(atoms))

    if isinstance(selection, Atom):
        return np.argwhere(np.array(atoms) == selection).flatten()

    if isinstance(selection, int):
        return np.array([selection])

    if isinstance(selection, slice):
        return np.arange(len(atoms))[selection]

    return np.array(selection)


def _traj_to_com_traj_by_residue(
    trajectory: Trajectory,
    selection,
    topology: Topology | None,
) -> Trajectory:
    """Center of mass trajectory of the residues of a topology."""

    topology = trajectory.topology if topology is None else topology

    if not topology.residues:
        logger.error(
            (
                "The topology has no residues. Use a topology with residues, "
                "e.g. from a restart file together with a moldescriptor file."
            ),
            exception=PQValueError,
        )

    first_frame = trajectory[0]

    if topology.n_atoms != first_frame.n_atoms:
        logger.error(
            (
                f"The topology has {topology.n_atoms} atoms, but the "
                f"trajectory has {first_frame.n_atoms} atoms per frame."
            ),
            exception=PQValueError,
        )

    indices = _selected_indices(first_frame.atoms, selection)
    group_ids = topology.residue_numbers[indices]

    # one atom per residue named by its selected atoms
    unique_ids, inverse = np.unique(group_ids, return_inverse=True)
    names = [[] for _ in unique_ids]
    for index, group_index in zip(indices, inverse.reshape(-1)):
        names[group_index].append(first_frame.atoms[index].name)
    com_atoms = [
        Atom(''.join(group_names), use_guess_element=False)
        for group_names in names
    ]

    com_traj = Trajectory()
    for frame in trajectory:
        if frame.n_atoms != first_frame.n_atoms:
            logger.error(
                "All frames need the same number of atoms for by_residue.",
                exception=PQValueError,
            )

        positions = center_of_mass_groups(
            frame.pos[indices],
            frame.atomic_masses[indices],
            group_ids,
            frame.cell,
        )
        com_traj.append(
            AtomicSystem(pos=positions, atoms=com_atoms, cell=frame.cell)
        )

    return com_traj
