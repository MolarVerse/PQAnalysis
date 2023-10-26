from PQAnalysis.selection.selection import Selection
from PQAnalysis.traj.trajectory import Trajectory


def traj2comtraj(trajectory, selection=None, group=None, map_name_to_element=None):

    if not isinstance(selection, Selection):
        selection = Selection(selection)

    com_traj = Trajectory()
    for frame in trajectory:
        frame = frame[selection]
        com_traj.append(frame.compute_com(group=None))

    return com_traj
