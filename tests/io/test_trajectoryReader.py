import pytest
import numpy as np

from . import pytestmark

from PQAnalysis.io import TrajectoryReader
from PQAnalysis.io.traj_file.exceptions import FrameReaderError
from PQAnalysis.traj import Frame
from PQAnalysis.core import Cell, Atom
from PQAnalysis.atomicSystem import AtomicSystem


class TestTrajectoryReader:
    @pytest.mark.usefixtures("tmpdir")
    def test__init__(self):
        with pytest.raises(FileNotFoundError) as exception:
            TrajectoryReader("tmp.xyz")
        assert str(exception.value) == "File tmp.xyz not found."

        open("tmp.xyz", "w")
        reader = TrajectoryReader("tmp.xyz")
        assert reader.filename == "tmp.xyz"
        assert reader.frames == []

    @pytest.mark.usefixtures("tmpdir")
    def test_read(self):

        file = open("tmp.xyz", "w")
        print("2 1.0 1.0 1.0", file=file)
        print("", file=file)
        print("h 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("h 1.0 0.0 0.0", file=file)
        print("o 0.0 1.0 1.0", file=file)
        file.close()

        reader = TrajectoryReader("tmp.xyz")

        traj = reader.read()

        cell = Cell(1.0, 1.0, 1.0)
        atoms = [Atom(atom) for atom in ["h", "o"]]

        print(traj[0].topology.atoms)
        print(traj[1].topology.atoms)
        print(traj[0].cell)
        print(traj[1].cell)

        frame1 = Frame(system=AtomicSystem(
            atoms=atoms, pos=np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]), cell=cell))
        frame2 = Frame(system=AtomicSystem(
            atoms=atoms, pos=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0]]), cell=cell))

        assert traj[0] == frame1
        # NOTE: here cell is not none because of the consecutive reading of frames
        # Cell will be taken from the previous frame
        assert traj[1] == frame2

        reader = TrajectoryReader("tmp.xyz", md_format="qmcfc")

        with pytest.raises(FrameReaderError) as exception:
            reader.read()
        assert str(
            exception.value) == "The first atom in one of the frames is not X. Please use pimd_qmcf (default) md engine instead"

        file = open("tmp.xyz", "w")
        print("2 1.0 1.0 1.0", file=file)
        print("", file=file)
        print("X 0.0 0.0 0.0", file=file)
        print("o 0.0 1.0 0.0", file=file)
        print("2", file=file)
        print("", file=file)
        print("X 1.0 0.0 0.0", file=file)
        print("o 0.0 1.0 1.0", file=file)
        file.close()

        reader = TrajectoryReader("tmp.xyz", md_format="qmcfc")

        traj = reader.read()

        print(traj[0].atoms[0].name.lower())
        print(traj[1].atoms[0].name)
        print(frame1.atoms[0].name)
        print(frame2.atoms[0].name)

        assert traj[0] == frame1[1]
        assert traj[1] == frame2[1]

        traj = reader.read()

        filenames = ["tmp.xyz", "tmp.xyz"]
        reader = TrajectoryReader(filenames, md_format="qmcfc")

        ref_traj = traj + traj
        traj = reader.read()

        assert traj == ref_traj
