import numpy as np

from PQAnalysis.io import RestartFileWriter
from PQAnalysis.traj import MDEngineFormat, Frame
from PQAnalysis.core import Cell, Atom
from PQAnalysis.atomicSystem import AtomicSystem
from PQAnalysis.topology import Topology


class TestRestartWriter:
    def test__init__(self):
        writer = RestartFileWriter("restart.dat")
        assert writer.filename == "restart.dat"
        assert writer.format == MDEngineFormat.PIMD_QMCF

    def test__write_box(self, capsys):
        writer = RestartFileWriter()
        cell = Cell(10.0, 10.0, 10.0)
        writer._write_box(cell)

        captured = capsys.readouterr()
        assert captured.out == "Box  10.0 10.0 10.0  90 90 90\n"

    def test__write_atoms(self, capsys):
        writer = RestartFileWriter()
        atoms = [Atom("C"), Atom("H"), Atom("H")]
        positions = np.array([[0.0, 0.0, 0.0],
                              [1.0, 1.0, 1.0],
                              [2.0, 2.0, 2.0]])
        frame = Frame(AtomicSystem(atoms, positions))

        print()
        writer._write_atoms(frame)

        captured = capsys.readouterr()
        assert captured.out == """
C    0    0    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
H    1    0    1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 
H    2    0    2.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 
"""

        frame = Frame(AtomicSystem(topology=Topology(
            atoms=atoms, residue_ids=np.array([1, 2, 3])), pos=positions))

        writer.format = MDEngineFormat.QMCFC

        print()
        writer._write_atoms(frame)

        captured = capsys.readouterr()
        assert captured.out == """
C    0    1    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
H    1    2    1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
H    2    3    2.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0
"""

    def test_write(self, capsys):
        writer = RestartFileWriter()
        atoms = [Atom("C"), Atom("H"), Atom("H")]
        positions = np.array([[0.0, 0.0, 0.0],
                              [1.0, 1.0, 1.0],
                              [2.0, 2.0, 2.0]])
        velocities = np.array([[0.0, 0.0, 0.0],
                              [1.0, 1.0, 1.0],
                              [2.0, 2.0, 2.0]])
        forces = np.array([[0.0, 0.0, 0.0],
                           [1.0, 1.0, 1.0],
                           [2.0, 2.0, 2.0]])
        cell = Cell(10.0, 10.0, 10.0)
        frame = Frame(AtomicSystem(atoms, positions, cell=cell,
                      vel=velocities, forces=forces))

        print()
        writer.write(frame)

        captured = capsys.readouterr()
        assert captured.out == """
Box  10.0 10.0 10.0  90 90 90
C    0    0    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
H    1    0    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 
H    2    0    2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 
"""
