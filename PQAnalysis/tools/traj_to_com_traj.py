"""
A module containing the tool to compute a center of mass trajectory for a given selection.
"""

from PQAnalysis.traj import Trajectory
from PQAnalysis.type_checking import runtime_type_checking



# TODO: add atom to element mapper if atom name not element names
# TODO rethink concept of selection/getitem/slices
@runtime_type_checking
def traj_to_com_traj(
    trajectory: Trajectory,
    selection=None,
    group=None
) -> Trajectory:
    """
    Function that computes the center of mass trajectory for a given selection.

    With the group parameter, the trajectory will be group according to the given group.
    For example if group = 2, the trajectory will be grouped in pairs of 2 and the center 
    of mass of each pair will be computed.

    Parameters
    ----------
    trajectory : Trajectory
        The trajectory to compute the center of mass trajectory for.
    selection : Selection, optional
        The selection to compute the center of mass trajectory for. 
        If None, the selection is all atoms.
    group : int, optional
        The group to compute the center of mass trajectory for.
        If None, the group is all atoms.
    """

    if len(trajectory) == 0:
        return Trajectory()

    com_traj = Trajectory()
    for frame in trajectory:

        if selection is None:
            selection = slice(0, frame.n_atoms)

        frame = frame[selection]
        com_traj.append(frame.compute_com_atomic_system(group=group))

    return com_traj
