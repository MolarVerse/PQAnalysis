import pytest

import numpy as np

from . import pytestmark

from PQAnalysis.io import RestartFileWriter
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.core import Cell, Atom
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.topology import Topology
from PQAnalysis.io.restart_file.exceptions import RestartFileWriterError



class TestRestartWriter:

    def test__init__(self):
        writer = RestartFileWriter("restart.dat")
        assert writer.filename == "restart.dat"
        assert writer.md_engine_format == MDEngineFormat.PQ

    def test__write_box(self, capsys):
        writer = RestartFileWriter()
        cell = Cell(10.0, 10.0, 10.0)
        lines = writer._get_box_line(cell)

        assert lines == "Box  10.0 10.0 10.0  90 90 90"

    def test__get_atom_lines(self, capsys):
        writer = RestartFileWriter()
        atoms = [Atom("C"), Atom("H"), Atom("H")]
        positions = np.array(
            [[0.0,
            0.0,
            0.0],
            [1.0,
            1.0,
            1.0],
            [2.0,
            2.0,
            2.0]]
        )

        frame = AtomicSystem(
            topology=Topology(atoms=atoms,
            residue_ids=np.array([1,
            2,
            3])),
            pos=positions
        )

        lines = writer._get_atom_lines(
            frame,
            md_engine_format=MDEngineFormat.QMCFC
        )

        assert lines == [
            "C    0    1    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0",
            "H    1    2    1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0",
            "H    2    3    2.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0",
        ]

        system = AtomicSystem(
            atoms=[Atom("C"),
            Atom("H"),
            Atom("H")],
            pos=np.array([[0.0,
            0.0,
            0.0],
            [1.0,
            1.0,
            1.0],
            [2.0,
            2.0,
            2.0]]),
        )

        atom_counter = np.array([0, 1])

        with pytest.raises(RestartFileWriterError) as exc:
            RestartFileWriter._get_atom_lines(system, atom_counter)
        assert str(exc.value) == (
            "The atom counter has to have the same length as "
            "the number of atoms in the frame if it is given as an array."
        )

        atom_lines = RestartFileWriter._get_atom_lines(system, 0)
        assert atom_lines[
            0] == "C    0    0    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
        assert atom_lines[
            1] == "H    0    0    1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0"
        assert atom_lines[
            2] == "H    0    0    2.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0"

        atom_lines = RestartFileWriter._get_atom_lines(
            system,
            np.array([0,
            1,
            2])
        )
        assert atom_lines[
            0] == "C    0    0    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
        assert atom_lines[
            1] == "H    1    0    1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0"
        assert atom_lines[
            2] == "H    2    0    2.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0"

    def test_write(self, capsys):
        writer = RestartFileWriter()
        atoms = [Atom("C"), Atom("H"), Atom("H")]
        positions = np.array(
            [[0.0,
            0.0,
            0.0],
            [1.0,
            1.0,
            1.0],
            [2.0,
            2.0,
            2.0]]
        )
        velocities = np.array(
            [[0.0,
            0.0,
            0.0],
            [1.0,
            1.0,
            1.0],
            [2.0,
            2.0,
            2.0]]
        )
        forces = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
        cell = Cell(10.0, 10.0, 10.0)
        frame = AtomicSystem(
            atoms,
            positions,
            cell=cell,
            vel=velocities,
            forces=forces
        )

        print()
        writer.write(frame)

        captured = capsys.readouterr()
        assert captured.out == """
Box  10.0 10.0 10.0  90 90 90
C    0    0    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
H    1    0    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
H    2    0    2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0
"""

        print()
        writer.md_engine_format = MDEngineFormat.QMCFC
        writer.write(frame)

        captured = capsys.readouterr()
        assert captured.out == """
Box  10.0 10.0 10.0  90 90 90
C    0    0    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
H    1    0    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
H    2    0    2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0
"""
