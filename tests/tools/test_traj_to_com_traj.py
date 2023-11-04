from PQAnalysis.tools.traj_to_com_traj import traj_to_com_traj

from PQAnalysis.traj.trajectory import Trajectory
from PQAnalysis.traj.frame import Frame


def test_traj_to_com_traj():
    traj = Trajectory()

    assert traj_to_com_traj(traj) == Trajectory()

    frame = Frame(atoms=['C', 'C', 'C'], coordinates=[
                  [0, 0, 0], [1, 1, 1], [2, 2, 2]])
    traj.append(frame)

    assert traj_to_com_traj(traj) == Trajectory(
        frames=[Frame(atoms=['ccc'], coordinates=[[1, 1, 1]])])

    frame = Frame(atoms=['C', 'C', 'C'], coordinates=[
                  [0, 0, 1], [1, 1, 2], [2, 2, 3]])

    traj.append(frame)

    assert traj_to_com_traj(traj) == Trajectory(frames=[Frame(
        atoms=['ccc'], coordinates=[[1, 1, 1]]), Frame(atoms=['ccc'], coordinates=[[1, 1, 2]])])

    assert traj_to_com_traj(traj, selection=[0, 1]) == Trajectory(frames=[Frame(
        atoms=['cc'], coordinates=[[0.5, 0.5, 0.5]]), Frame(atoms=['cc'], coordinates=[[0.5, 0.5, 1.5]])])

    traj = Trajectory()
    frame = Frame(atoms=['C', 'C', 'H', 'H'], coordinates=[
                  [0, 0, 1], [1, 1, 2], [2, 2, 3], [3, 3, 4]])
    traj.append(frame)
    frame = Frame(atoms=['C', 'C', 'H', 'H'], coordinates=[
                  [0, 1, 1], [1, 2, 2], [2, 3, 3], [3, 4, 4]])
    traj.append(frame)

    assert traj_to_com_traj(traj, group=2) == Trajectory(frames=[Frame(atoms=['cc', 'hh'], coordinates=[
        [0.5, 0.5, 1.5], [2.5, 2.5, 3.5]]), Frame(atoms=['cc', 'hh'], coordinates=[[0.5, 1.5, 1.5], [2.5, 3.5, 3.5]])])
