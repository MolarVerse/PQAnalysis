from PQAnalysis.selection.selection import Selection
from PQAnalysis.traj.trajectory import Trajectory


def traj2comtraj(selection, trajectory, group=None, map_name_to_element=None):

    if not isinstance(selection, Selection):
        selection = Selection(selection)

    if group != None:
        if len(selection) % group != 0:
            raise ValueError('Invalid group size.')
        else:
            selection = selection.reshape(len(selection) // group, group)
    else:
        selection = selection.reshape(len(selection), 1)

    com_traj = Trajectory()
    for frame in trajectory:
        com_traj.append(frame.compute_com(selection))

    return com_traj
