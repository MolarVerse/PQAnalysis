"""
A module including the API for the topology subpackage.
"""

from beartype.typing import List

from PQAnalysis.traj import Trajectory
from PQAnalysis.io import FileWritingMode, read_restart_file
from PQAnalysis.types import Np1DIntArray
from PQAnalysis.type_checking import runtime_type_checking

from .selection import SelectionCompatible, Selection
from .shake_topology import ShakeTopologyGenerator



@runtime_type_checking
def generate_shake_topology_file(
    trajectory: Trajectory,
    selection: SelectionCompatible = None,
    output: str | None = None,
    mode: FileWritingMode | str = "w",
    use_full_atom_info: bool = False,
) -> None:
    """
    Wrapper function to generate a shake topology file for a given trajectory.

    Parameters
    ----------
    trajectory : Trajectory
        The trajectory to generate the shake topology for.
    selection : SelectionCompatible, optional
        Selection is either a selection object or any object that can be
        initialized via 'Selection(selection)'. default None (all atoms)
    output : str | None, optional
        The output file. If not specified, the output is printed to stdout.
    mode : FileWritingMode | str, optional
        The writing mode, by default "w". The following modes are available:
        - "w": write
        - "a": append
        - "o": overwrite
    use_full_atom_info : bool, optional
        If True, the full atom information (name, index, mass) is used
        for the selection, by default False
        Is always ignored if atoms is not a list of atom objects.
    """

    generator = ShakeTopologyGenerator(selection, use_full_atom_info)
    generator.generate_topology(trajectory)
    generator.write_topology(output, mode)



@runtime_type_checking
def select_from_restart_file(
    selection: SelectionCompatible,
    restart_file: str,
    moldescriptor_file: str | None = None,
    use_full_atom_info: bool = False,
) -> Np1DIntArray:
    """
    Selects atoms from a restart file and writes them to a new file.

    Parameters
    ----------
    selection : SelectionCompatible
        Selection is either a selection object or any object that can be
        initialized via 'Selection(selection)'.
    restart_file : str
        The restart file to read the atoms from.
    moldescriptor_file : str | None, optional
        The moldescriptor file to read the atom types from, by default None
    use_full_atom_info : bool, optional
        If True, the full atom information (name, index, mass) is used
        for the selection, by default False
        Is always ignored if atoms is not a list of atom objects.
    """

    system = read_restart_file(restart_file, moldescriptor_file)

    selection = Selection(selection)

    return selection.select(
        system.topology,
        use_full_atom_info=use_full_atom_info
    )



def selection_from_restart_file_as_list(
    selection: SelectionCompatible,
    restart_file: str,
    moldescriptor_file: str | None = None,
    use_full_atom_info: bool = False,
) -> List[int]:
    """
    Selects atoms from a restart file and writes them to a new file.

    Parameters
    ----------
    selection : SelectionCompatible
        Selection is either a selection object or any object that can be
        initialized via 'Selection(selection)'.
    restart_file : str
        The restart file to read the atoms from.
    moldescriptor_file : str | None, optional
        The moldescriptor file to read the atom types from, by default None
    use_full_atom_info : bool, optional
        If True, the full atom information (name, index, mass) is used
        for the selection, by default False
        Is always ignored if atoms is not a list of atom objects.
    """

    return list(
        select_from_restart_file(
        selection,
        restart_file,
        moldescriptor_file,
        use_full_atom_info
        )
    )
