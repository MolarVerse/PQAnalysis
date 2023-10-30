"""
A module containing the tool to compute a center of mass trajectory for a given selection.
"""

from typing import Union, List

from PQAnalysis.traj.selection import Selection
from PQAnalysis.traj.trajectory import Trajectory


# TODO: add atom to element mapper if atomname not element names
def traj_to_com_traj(trajectory: Trajectory, selection: Union[List[int], Selection] = None, group=None):
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
        The selection to compute the center of mass trajectory for. If None, the selection is all atoms.
    group : int, optional
        The group to compute the center of mass trajectory for. If None, the group is all atoms.
    """
    if not isinstance(selection, Selection):
        selection = Selection(selection)

    com_traj = Trajectory()
    for frame in trajectory:
        frame = frame[selection]
        com_traj.append(frame.compute_com(group=group))

    return com_traj
