from .selection import SelectionCompatible
from .shakeTopology import ShakeTopologyGenerator
from PQAnalysis.traj import Trajectory
from PQAnalysis.io import FileWritingMode


def generate_shake_topology_file(trajectory: Trajectory,
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
        Selection is either a selection object or any object that can be initialized via 'Selection(selection)'. default None (all atoms)
    output : str | None, optional
        The output file. If not specified, the output is printed to stdout.
    mode : FileWritingMode | str, optional
        The writing mode, by default "w". The following modes are available:
        - "w": write
        - "a": append
        - "o": overwrite
    use_full_atom_info : bool, optional
        If True, the full atom information (name, index, mass) is used for the selection, by default False
        Is always ignored if atoms is not a list of atom objects.
    """

    generator = ShakeTopologyGenerator(selection, use_full_atom_info)
    generator.generate_topology(trajectory)
    generator.write_topology(output, mode)
