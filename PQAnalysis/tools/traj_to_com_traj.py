"""
A module containing the tool to compute a center of mass trajectory for a given selection.
"""

import logging

import numpy as np

from beartype.typing import List

from PQAnalysis import __package_name__
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.atomic_system._center_of_mass import center_of_mass_groups
from PQAnalysis.core import Atom
from PQAnalysis.exceptions import PQValueError
from PQAnalysis.io import TrajectoryReader, TrajectoryWriter, read_restart_file
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.topology import Selection, SelectionCompatible, Topology
from PQAnalysis.traj import MDEngineFormat, Trajectory
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.types import PositiveInt
from PQAnalysis.utils.custom_logging import setup_logger

logger = logging.getLogger(__package_name__).getChild("traj_to_com_traj")
logger = setup_logger(logger)

#: Elements from this atomic number on (X, Q, Cav, Sup) are placeholders
#: with a made-up mass, not real atoms.
PLACEHOLDER_ATOMIC_NUMBER = 999



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
    by the residue number and named by the residue name of the molecule
    (e.g. the name given in the moldescriptor file). It carries no residues
    itself. A molecule has to be smaller than half of the box in every
    direction; an error is raised for molecules whose atoms spread over half
    of the box or more around the first atom (larger molecules in which this
    is not visible, such as a pair of atoms farther apart than half the box,
    cannot be detected). A warning is logged if selected
    atoms have a placeholder mass (elements X, Q, Cav, Sup), because they
    dominate the center of mass of their molecule.

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
    AtomicSystemError
        With by_residue, if a molecule spans half of the box or more.
    AtomicSystemMassError
        If a selected atom has no mass.
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

    _warn_about_placeholder_atoms(first_frame.atoms, indices)

    # one atom per residue, named by the residue name
    com_atoms = [
        Atom(topology.residues[residue_number].name, use_guess_element=False)
        for residue_number in np.unique(group_ids)
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
            check_extent=True,
        )
        com_traj.append(
            AtomicSystem(pos=positions, atoms=com_atoms, cell=frame.cell)
        )

    return com_traj


def _warn_about_placeholder_atoms(atoms, indices) -> None:
    """Warns if selected atoms have a placeholder mass (X, Q, Cav, Sup)."""

    placeholders = sorted({
        atoms[index].name
        for index in indices
        if (atoms[index].atomic_number or 0) >= PLACEHOLDER_ATOMIC_NUMBER
    })

    if placeholders:
        logger.warning(
            (
                f"The selection contains atoms with a placeholder mass "
                f"({', '.join(placeholders)}). They dominate the center of "
                "mass of their molecule; exclude them with the selection."
            )
        )



@runtime_type_checking
def traj2com(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    trajectory_files: List[str],
    output: str | None = None,
    mode: FileWritingMode | str = "w",
    md_format: MDEngineFormat | str = MDEngineFormat.PQ,
    selection: SelectionCompatible = None,
    use_full_atom_info: bool = False,
    group: PositiveInt | None = None,
    by_residue: bool = False,
    restart_file: str | None = None,
    moldescriptor_file: str | None = None,
) -> None:
    """
    Reads xyz trajectory files, computes the center of mass trajectory with
    :py:func:`traj_to_com_traj` and writes it to a file or to stdout.

    For by_residue, the molecules are defined by a restart file together with
    the moldescriptor file, as for the other analyses: the restart file gives
    the molecule of every atom and the moldescriptor file the molecule types.
    The molecule types may contain the WATER_TYPE and AMMONIA_TYPE lines of PQ.

    Parameters
    ----------
    trajectory_files : list of str
        The xyz trajectory file(s) to read.
    output : str, optional
        The output file. If not specified, the output is printed to stdout.
    mode : FileWritingMode | str, optional
        The writing mode, by default "w". The following modes are available:
        - "w": write
        - "a": append
        - "o": overwrite
    md_format : MDEngineFormat | str, optional
        The format of the input trajectory. The default is MDEngineFormat.PQ.
    selection : SelectionCompatible, optional
        The atoms that define the centers of mass, as a selection string
        (see :py:class:`~PQAnalysis.topology.selection.Selection`). If None,
        all atoms are selected.
    use_full_atom_info : bool, optional
        Whether the selection uses the full atom information (atom type
        names and elements) or the elements only, by default False.
    group : PositiveInt, optional
        The number of consecutive atoms for each center of mass. If neither
        group nor by_residue is given, the center of mass of all selected
        atoms is computed.
    by_residue : bool, optional
        Compute the center of mass of every molecule, by default False. Needs
        restart_file and moldescriptor_file. Cannot be combined with group.
    restart_file : str, optional
        The restart file that defines the molecules of the system.
    moldescriptor_file : str, optional
        The moldescriptor file with the molecule types, needs restart_file.

    Raises
    ------
    PQValueError
        If moldescriptor_file is given without restart_file, if the selection
        contains no atoms or for the errors of :py:func:`traj_to_com_traj`.
    """

    if moldescriptor_file is not None and restart_file is None:
        logger.error(
            "The moldescriptor file needs a restart file that defines the "
            "molecules.",
            exception=PQValueError,
        )

    topology = None
    if restart_file is not None:
        topology = read_restart_file(
            restart_file,
            moldescriptor_file,
            md_engine_format=md_format,
        ).topology

    trajectory = TrajectoryReader(
        trajectory_files,
        md_format=md_format,
        topology=topology,
    ).read()

    indices = None
    if selection is not None and len(trajectory) > 0:
        indices = Selection(selection).select(
            trajectory.topology, use_full_atom_info
        )

        if len(indices) == 0:
            logger.error(
                "The selection does not select any atoms.",
                exception=PQValueError,
            )

    com_traj = traj_to_com_traj(
        trajectory,
        selection=indices,
        group=group,
        by_residue=by_residue,
    )

    TrajectoryWriter(
        filename=output,
        engine_format=md_format,
        mode=mode,
    ).write(com_traj)
