import numpy as np

from PQAnalysis.tools.traj_to_com_traj import traj_to_com_traj

from PQAnalysis.traj.trajectory import Trajectory
from PQAnalysis.traj.frame import Frame
from PQAnalysis.core.atom import Atom
from PQAnalysis.core.atomicSystem import AtomicSystem


def test_traj_to_com_traj():
    traj = Trajectory()

    assert traj_to_com_traj(traj) == Trajectory()

    atoms = [Atom(atom) for atom in ['C', 'C', 'C']]
    coordinates = [[0, 0, 0], [1, 1, 1], [2, 2, 2]]

    frame = Frame(AtomicSystem(atoms=atoms, pos=np.array(coordinates)))
    traj.append(frame)

    print(traj[0].n_atoms)

    traj_output = traj_to_com_traj(traj)

    assert np.allclose(traj_output[0].pos, [[1, 1, 1]])
    assert traj_output[0].system.combined_name == 'CCC'

    coordinates = [[0, 0, 1], [1, 1, 2], [2, 2, 3]]
    frame = Frame(AtomicSystem(atoms=atoms, pos=np.array(coordinates)))

    traj.append(frame)

    traj_output = traj_to_com_traj(traj)

    assert np.allclose(traj_output[0].pos, [[1, 1, 1]])
    assert traj_output[0].system.combined_name == 'CCC'
    assert np.allclose(traj_output[1].pos, [[1, 1, 2]])
    assert traj_output[1].system.combined_name == 'CCC'

    traj_output = traj_to_com_traj(traj, selection=slice(0, 2))
    assert np.allclose(traj_output[0].pos, [[0.5, 0.5, 0.5]])
    assert traj_output[0].system.combined_name == 'CC'
    assert np.allclose(traj_output[1].pos, [[0.5, 0.5, 1.5]])
    assert traj_output[1].system.combined_name == 'CC'

    traj = Trajectory()
    atoms = [Atom(atom) for atom in ['C', 'C', 'H', 'H']]
    coordinates1 = [[0, 0, 1], [1, 1, 2], [2, 2, 3], [3, 3, 4]]
    coordinates2 = [[0, 1, 1], [1, 2, 2], [2, 3, 3], [3, 4, 4]]
    frame1 = Frame(AtomicSystem(atoms=atoms, pos=np.array(coordinates1)))
    frame2 = Frame(AtomicSystem(atoms=atoms, pos=np.array(coordinates2)))
    traj.append(frame1)
    traj.append(frame2)

    traj_output = traj_to_com_traj(traj, group=2)

    assert np.allclose(traj_output[0].pos, [[0.5, 0.5, 1.5], [2.5, 2.5, 3.5]])
    assert traj_output[0].atoms == [
        Atom('CC', use_guess_element=False), Atom('HH', use_guess_element=False)]
    assert traj_output[0].system.combined_name == 'CCHH'
    assert np.allclose(traj_output[1].pos, [[0.5, 1.5, 1.5], [2.5, 3.5, 3.5]])
    assert traj_output[1].atoms == [
        Atom('CC', use_guess_element=False), Atom('HH', use_guess_element=False)]
    assert traj_output[1].system.combined_name == 'CCHH'
