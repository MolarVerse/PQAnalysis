"""
A module containing the tool to compute a center of mass trajectory for a given selection.
"""

import numpy as np

from beartype.typing import List

from PQAnalysis.traj.trajectory import Trajectory


# TODO: add atom to element mapper if atomname not element names
# TODO rethink concept of selection/getitem/slices
def traj_to_com_traj(trajectory: Trajectory, selection=None, group=None) -> Trajectory:
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

    if len(trajectory) == 0:
        return Trajectory()

    com_traj = Trajectory()
    for frame in trajectory:

        if selection is None:
            selection = slice(0, frame.n_atoms)

        frame = frame[selection]
        com_traj.append(frame.compute_com_frame(group=group))

    return com_traj
