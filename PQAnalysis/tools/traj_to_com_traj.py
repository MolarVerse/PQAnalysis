"""
A module containing the tool to compute a center of mass trajectory for a given selection.
"""

import logging

import numpy as np

from beartype.typing import List, Sequence

from PQAnalysis import __package_name__
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.atomic_system._center_of_mass import (
    center_of_mass_properties,
    placeholder_atom_names,
)
from PQAnalysis.core import Atom
from PQAnalysis.exceptions import PQValueError
from PQAnalysis.io import (
    MoldescriptorReader,
    TrajectoryReader,
    TrajectoryWriter,
    read_restart_file,
)
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.topology import Selection, SelectionCompatible, Topology
from PQAnalysis.traj import MDEngineFormat, Trajectory
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.types import PositiveInt
from PQAnalysis.utils.custom_logging import setup_logger

logger = logging.getLogger(__package_name__).getChild("traj_to_com_traj")
logger = setup_logger(logger)



# TODO: add atom to element mapper if atom name not element names
# TODO rethink concept of selection/getitem/slices
@runtime_type_checking
def traj_to_com_traj(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    trajectory: Trajectory,
    selection=None,
    group=None,
    by_residue: bool = False,
    topology: Topology | None = None,
    continuous: bool = False,
    molecule_types: Sequence[str | int] | None = None,
) -> Trajectory:
    """
    Function that computes the center of mass trajectory for a given selection.

    The velocities, forces and charges of the frames, if they have them, are
    combined as well: a center of mass has the mass weighted mean of the
    velocities of its atoms, the sum of their forces and the sum of their
    charges. A trajectory that only has velocities (a vel file) gives a
    trajectory of center of mass velocities, which can be used for the
    velocity autocorrelation function.

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
    continuous : bool, optional
        If True, the positions of the centers of mass follow their way over
        the periodic boundary from frame to frame instead of being folded
        back into the cell, so that they can be used for the mean square
        displacement. Every frame is placed at the image of the center that
        is closest to the previous frame (in fractional coordinates, so a
        cell that changes in size is handled), so a center must not move
        farther than half of the box between two frames. The first frame is
        folded into the cell. Has no effect in a vacuum cell. By default
        False.
    molecule_types : Sequence[str | int], optional
        With by_residue, only the molecules of these types get a center of
        mass, e.g. ``["TIP3P-Water"]`` or ``[1]``. A type is given by the name
        of the mol type in the moldescriptor file or by its number (the id of
        the residue). Combined with the selection, the center of a molecule
        is computed from the selected atoms of the molecules of these types.
        By default None (all molecules).

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
    PQValueError
        If continuous is used for a trajectory without positions, if
        molecule_types is used without by_residue, or if a molecule type
        does not exist in the topology.
    """

    if by_residue and group is not None:
        logger.error(
            "The parameters group and by_residue cannot be combined.",
            exception=PQValueError,
        )

    if molecule_types is not None and not by_residue:
        logger.error(
            "The parameter molecule_types needs by_residue.",
            exception=PQValueError,
        )

    if len(trajectory) == 0:
        return Trajectory()

    if continuous and not trajectory[0].has_pos:
        logger.error(
            "A continuous center of mass needs the positions.",
            exception=PQValueError,
        )

    if by_residue:
        com_traj = _traj_to_com_traj_by_residue(
            trajectory, selection, topology, molecule_types
        )
    else:
        com_traj = Trajectory()
        for frame in trajectory:
            selected = frame if selection is None else frame[selection]
            com_traj.append(selected.compute_com_atomic_system(group=group))

    if continuous:
        _make_continuous(com_traj)

    return com_traj


def _make_continuous(com_traj: Trajectory) -> None:
    """
    Places every frame at the periodic image closest to the previous frame.

    The first frame stays as it is. The images are followed in fractional
    coordinates, so that a cell that changes in size from frame to frame is
    handled correctly. A vacuum cell has no periodic images, nothing changes.
    """

    previous = None

    for frame in com_traj:
        cell = frame.cell

        if cell.is_vacuum:
            previous = None
            continue

        fractional = frame.pos @ cell.inverse_box_matrix.T

        if previous is not None:
            step = fractional - previous
            fractional = previous + step - np.round(step)
            frame.pos = fractional @ cell.box_matrix.T

        previous = fractional


def _in_molecule_types(
    topology: Topology, molecule_types: Sequence[str | int]
) -> np.ndarray:
    """
    Returns for every residue of the topology whether it is of one of the
    given types, which are given by name or by number.
    """

    names = {m for m in molecule_types if isinstance(m, str)}
    numbers = {m for m in molecule_types if not isinstance(m, str)}
    available_names = sorted({residue.name for residue in topology.residues})
    available_numbers = sorted({residue.id for residue in topology.residues})

    unknown = [
        m for m in molecule_types
        if (m not in available_names and m not in available_numbers)
    ]
    if unknown:
        logger.error(
            (
                f"The molecule types {unknown} are not in the topology, which "
                f"has the types {available_names} with the numbers "
                f"{available_numbers}."
            ),
            exception=PQValueError,
        )

    return np.array([
        residue.name in names or residue.id in numbers
        for residue in topology.residues
    ])


def _take(values: np.ndarray, indices: np.ndarray) -> np.ndarray | None:
    """Returns the values of the given atoms, or None if there are none."""

    return values[indices] if len(values) > 0 else None


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
    molecule_types: Sequence[str | int] | None = None,
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

    if molecule_types is not None:
        indices = indices[
            _in_molecule_types(topology, molecule_types)[
                topology.residue_numbers[indices]]
        ]

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

        properties = center_of_mass_properties(
            frame.atomic_masses[indices],
            group_ids,
            frame.cell,
            positions=_take(frame.pos, indices),
            velocities=_take(frame.vel, indices),
            forces=_take(frame.forces, indices),
            charges=_take(frame.charges, indices),
            check_extent=True,
        )
        com_traj.append(
            AtomicSystem(atoms=com_atoms, cell=frame.cell, **properties)
        )

    return com_traj


def _warn_about_placeholder_atoms(atoms, indices) -> None:
    """Warns if selected atoms have a placeholder mass (X, Q, Cav, Sup)."""

    placeholders = placeholder_atom_names(atoms, indices)

    if placeholders:
        logger.warning(
            (
                f"The selection contains atoms with a placeholder mass "
                f"({', '.join(placeholders)}). They dominate the center of "
                "mass of their molecule; exclude them with the selection."
            )
        )



@runtime_type_checking
# pylint: disable-next=too-many-locals
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
    continuous: bool = False,
    molecule_types: Sequence[str | int] | None = None,
) -> None:
    """
    Reads trajectory files, computes the center of mass trajectory with
    :py:func:`traj_to_com_traj` and writes it to a file or to stdout.

    The format of the trajectory is inferred from the file extension (xyz,
    vel, force or charge) and the result is written in the same format. A vel
    trajectory therefore gives the velocities of the centers of mass, a force
    trajectory their forces and a charge trajectory their charges.

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
    continuous : bool, optional
        Follow the centers over the periodic boundary instead of folding them
        back into the cell, see :py:func:`traj_to_com_traj`. By default False.
    molecule_types : Sequence[str | int], optional
        With by_residue, only the molecules of these types get a center of
        mass, see :py:func:`traj_to_com_traj`. In addition to the names and
        numbers of the mol types, the keywords ``water`` and ``ammonia``
        (any case) select the mol types that the moldescriptor file names
        with its WATER_TYPE and AMMONIA_TYPE lines.

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

    molecule_types = _resolve_type_keywords(molecule_types, moldescriptor_file)

    topology = None
    if restart_file is not None:
        topology = read_restart_file(
            restart_file,
            moldescriptor_file,
            md_engine_format=md_format,
        ).topology

    reader = TrajectoryReader(
        trajectory_files,
        md_format=md_format,
        topology=topology,
    )
    trajectory = reader.read()

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
        continuous=continuous,
        molecule_types=molecule_types,
    )

    TrajectoryWriter(
        filename=output,
        engine_format=md_format,
        mode=mode,
    ).write(com_traj, traj_type=reader.traj_format)



def _resolve_type_keywords(
    molecule_types: Sequence[str | int] | None,
    moldescriptor_file: str | None,
) -> Sequence[str | int] | None:
    """
    Replaces the keywords water and ammonia by the mol type numbers that the
    moldescriptor file names with its WATER_TYPE and AMMONIA_TYPE lines.
    """

    if molecule_types is None:
        return None

    keywords = {"water": "water_type", "ammonia": "ammonia_type"}
    resolved = []
    reader = None

    for molecule_type in molecule_types:
        attribute = keywords.get(str(molecule_type).lower())

        if attribute is None:
            resolved.append(molecule_type)
            continue

        if moldescriptor_file is None:
            logger.error(
                (
                    f"The molecule type '{molecule_type}' is a keyword of the "
                    "WATER_TYPE and AMMONIA_TYPE lines and needs the "
                    "moldescriptor file."
                ),
                exception=PQValueError,
            )

        if reader is None:
            reader = MoldescriptorReader(moldescriptor_file)
            reader.read()

        number = getattr(reader, attribute)
        if number is None:
            logger.error(
                (
                    f"The moldescriptor file '{moldescriptor_file}' has no "
                    f"{attribute.upper()} line for the molecule type "
                    f"'{molecule_type}'."
                ),
                exception=PQValueError,
            )

        resolved.append(number)

    return resolved
