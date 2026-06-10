import pytest
import sys
import numpy as np

from . import pytestmark

from PQAnalysis.io import (
    ExtXYZProfile,
    TrajectoryWriter,
    read_trajectory,
    write_trajectory,
    FileWritingMode,
)
from PQAnalysis.traj import Trajectory, TrajectoryFormat, MDEngineFormat
from PQAnalysis.core import Cell, Atom
from PQAnalysis.io.exceptions import ExtXYZProfileError
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.traj.exceptions import MDEngineFormatError
from PQAnalysis.utils.units import eV, kcal_per_mol



def test_write_trajectory(capsys):
    atoms = [Atom(atom) for atom in ['h', 'o']]
    coordinates1 = np.array([[0, 0, 0], [0, 0, 1]])
    coordinates2 = np.array([[0, 0, 0], [0, 0, 1]])
    frame1 = AtomicSystem(atoms=atoms, pos=coordinates1, cell=Cell(10, 10, 10))
    frame2 = AtomicSystem(atoms=atoms, pos=coordinates2, cell=Cell(11, 10, 10))
    traj = Trajectory([frame1, frame2])

    print()
    write_trajectory(traj, engine_format="PQ")

    captured = capsys.readouterr()
    assert captured.out == """
2 10 10 10 90 90 90

h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
2 11 10 10 90 90 90

h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
"""


@pytest.mark.usefixtures("tmpdir")
def test_write_extxyz_roundtrip():
    atoms = [Atom("h"), Atom("o")]
    virial = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
    stress = np.array([[9.0, 8.0, 7.0], [6.0, 5.0, 4.0], [3.0, 2.0, 1.0]])
    frame = AtomicSystem(
        atoms=atoms,
        pos=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
        vel=np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]),
        forces=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        charges=np.array([-0.1, -0.2]),
        energy=-1.5,
        virial=virial,
        stress=stress,
        cell=Cell(10.0, 11.0, 12.0, 80.0, 85.0, 95.0),
    )

    write_trajectory(
        Trajectory([frame]),
        filename="output.extxyz",
        traj_type="extxyz",
    )

    with open("output.extxyz", "r", encoding="utf-8") as file:
        lines = file.read().splitlines()

    assert lines[0] == "2"
    assert "Lattice=" in lines[1]
    assert (
        'Properties="species:S:1:pos:R:3:vel:R:3:forces:R:3:charge:R:1"'
        in lines[1]
    )
    assert "energy=-1.5000000000" in lines[1]
    assert 'virial="1.0000000000 4.0000000000 7.0000000000' in lines[1]
    assert 'stress="9.0000000000 6.0000000000 3.0000000000' in lines[1]

    trajectory = read_trajectory("output.extxyz", traj_format="extxyz")
    output_frame = trajectory[0]

    assert output_frame.atoms == atoms
    assert np.allclose(output_frame.pos, frame.pos)
    assert np.allclose(output_frame.vel, frame.vel)
    assert np.allclose(output_frame.forces, frame.forces)
    assert np.allclose(output_frame.charges, frame.charges)
    assert output_frame.energy == frame.energy
    assert np.allclose(output_frame.virial, frame.virial)
    assert np.allclose(output_frame.stress, frame.stress)
    assert output_frame.cell == frame.cell


@pytest.mark.usefixtures("tmpdir")
def test_write_extxyz_ase_profile_converts_energy_like_units():
    from ase.io import read as ase_read

    atoms = [Atom("h"), Atom("o")]
    frame = AtomicSystem(
        atoms=atoms,
        pos=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
        forces=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        energy=-1.5,
        virial=np.array(
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
        ),
        stress=np.array(
            [[9.0, 8.0, 7.0], [6.0, 5.0, 4.0], [3.0, 2.0, 1.0]]
        ),
        cell=Cell(10.0, 11.0, 12.0),
    )
    conversion_factor = eV / kcal_per_mol

    writer = TrajectoryWriter("ase.extxyz", extxyz_profile=ExtXYZProfile.ASE)
    writer.write(frame, traj_type="extxyz")

    with open("ase.extxyz", "r", encoding="utf-8") as file:
        lines = file.read().splitlines()

    assert f"energy={-1.5 * conversion_factor:.10f}" in lines[1]
    assert f"{1.0 * conversion_factor:.10f}" in lines[2]
    assert f"{9.0 * conversion_factor:.10f}" in lines[1]

    ase_frame = ase_read("ase.extxyz", index=0)
    assert np.isclose(ase_frame.calc.results["energy"], -1.5 * conversion_factor)
    assert np.allclose(
        ase_frame.calc.results["forces"],
        frame.forces * conversion_factor,
    )


@pytest.mark.usefixtures("tmpdir")
def test_write_extxyz_nep_profile_uses_gpumd_convention():
    atoms = [Atom("h"), Atom("o")]
    frame = AtomicSystem(
        atoms=atoms,
        pos=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
        forces=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        energy=-1.5,
        virial=np.array(
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
        ),
        stress=np.array(
            [[9.0, 8.0, 7.0], [6.0, 5.0, 4.0], [3.0, 2.0, 1.0]]
        ),
        cell=Cell(10.0, 11.0, 12.0),
    )
    conversion_factor = eV / kcal_per_mol

    writer = TrajectoryWriter("nep.extxyz", extxyz_profile=ExtXYZProfile.NEP)
    writer.write(frame, traj_type="extxyz")

    with open("nep.extxyz", "r", encoding="utf-8") as file:
        lines = file.read().splitlines()

    assert lines[1].startswith('lattice="10.0000000000 0.0000000000')
    assert 'properties="species:S:1:pos:R:3:forces:R:3"' in lines[1]
    assert "Lattice=" not in lines[1]
    assert "Properties=" not in lines[1]
    assert f"energy={-1.5 * conversion_factor:.10f}" in lines[1]
    assert (
        f'virial="{1.0 * conversion_factor:.10f} '
        f"{2.0 * conversion_factor:.10f} "
        f"{3.0 * conversion_factor:.10f}"
    ) in lines[1]
    assert (
        f'stress="{9.0 * conversion_factor:.10f} '
        f"{8.0 * conversion_factor:.10f} "
        f"{7.0 * conversion_factor:.10f}"
    ) in lines[1]
    assert (
        f"H 0.0000000000 0.0000000000 0.0000000000 "
        f"{1.0 * conversion_factor:.10f} "
        f"{2.0 * conversion_factor:.10f} "
        f"{3.0 * conversion_factor:.10f}"
    ) == lines[2]



class TestTrajectoryWriter:

    def test__init__(self):

        with pytest.raises(MDEngineFormatError) as exception:
            TrajectoryWriter(engine_format="notAFormat")
        assert str(exception.value) == (
            "\n"
            "'notaformat' is not a valid MDEngineFormat.\n"
            f"Possible values are: {MDEngineFormat.member_repr()} "
            "or their case insensitive string representation: "
            f"{MDEngineFormat.value_repr()}"
        )

        writer = TrajectoryWriter()
        assert writer.file == sys.stdout
        assert writer.filename is None
        assert writer.mode == FileWritingMode.WRITE
        assert writer.format == MDEngineFormat.PQ
        assert writer.extxyz_profile == ExtXYZProfile.PQ

        writer = TrajectoryWriter(engine_format="qmcfc")
        assert writer.format == MDEngineFormat.QMCFC

        writer = TrajectoryWriter(engine_format="PQ")
        assert writer.format == MDEngineFormat.PQ

        writer.extxyz_profile = "ASE"
        assert writer.extxyz_profile == ExtXYZProfile.ASE

        with pytest.raises(ExtXYZProfileError) as exception:
            TrajectoryWriter(extxyz_profile="notAProfile")
        assert str(exception.value) == (
            "\n"
            "'notaprofile' is not a valid ExtXYZProfile.\n"
            f"Possible values are: {ExtXYZProfile.member_repr()} "
            "or their case insensitive string representation: "
            f"{ExtXYZProfile.value_repr()}"
        )

    def test__write_header(self, capsys):

        writer = TrajectoryWriter()
        assert writer.mode == FileWritingMode.WRITE
        writer._write_header(1, Cell(10, 10, 10))

        captured = capsys.readouterr()
        assert captured.out == "1 10 10 10 90 90 90\n"

        writer._write_header(1)
        captured = capsys.readouterr()
        assert captured.out == "1\n"

        writer.format = MDEngineFormat.QMCFC
        writer._write_header(1)
        captured = capsys.readouterr()
        assert captured.out == "2\n"

    def test__write_comment(self, capsys):

        writer = TrajectoryWriter()
        writer.type = TrajectoryFormat.XYZ
        writer._write_comment(
            AtomicSystem(
            atoms=[Atom(atom) for atom in ["h",
            "o"]],
            cell=Cell(10,
            10,
            10)
            )
        )

        captured = capsys.readouterr()
        assert captured.out == "\n"

        forces = np.array([[1, 0, 3], [0, 2, 1]])
        writer.type = TrajectoryFormat.FORCE
        writer._write_comment(
            AtomicSystem(
            atoms=[Atom(atom) for atom in ["h",
            "o"]],
            cell=Cell(10,
            10,
            10),
            forces=forces
            )
        )

        captured = capsys.readouterr()
        assert captured.out == "sum of forces: 1.000000e+00 2.000000e+00 4.000000e+00\n"

    def test__write_xyz(self, capsys):

        writer = TrajectoryWriter()
        writer.type = TrajectoryFormat.XYZ

        print()
        writer._write_xyz(
            atoms=[Atom(atom) for atom in ["h",
            "o"]],
            xyz=np.array([[0,
            0,
            0],
            [0,
            0,
            1]])
        )

        captured = capsys.readouterr()
        assert captured.out == """
h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
"""

        writer.format = "qmcfc"

        print()
        writer._write_xyz(
            atoms=[Atom(atom) for atom in ["h",
            "o"]],
            xyz=np.array([[0,
            0,
            0],
            [0,
            0,
            1]])
        )

        captured = capsys.readouterr()
        assert captured.out == """
X   0.0 0.0 0.0
h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
"""

    def test__write_scalar(self, capsys):

        writer = TrajectoryWriter()
        writer._write_scalar(
            atoms=[Atom(atom) for atom in ["h",
            "o"]],
            scalar=np.array([1,
            2])
        )

        captured = capsys.readouterr()
        assert captured.out == "h 1\no 2\n"

    def test_write(self, capsys):

        atoms = [Atom(atom) for atom in ['h', 'o']]
        coordinates1 = np.array([[0, 0, 0], [0, 0, 1]])
        coordinates2 = np.array([[0, 0, 0], [0, 0, 1]])

        frame1 = AtomicSystem(
            atoms=atoms,
            pos=coordinates1,
            cell=Cell(10,
            10,
            10)
        )
        frame2 = AtomicSystem(
            atoms=atoms,
            pos=coordinates2,
            cell=Cell(11,
            10,
            10)
        )

        traj = Trajectory([frame1, frame2])
        writer = TrajectoryWriter()
        assert writer.mode == FileWritingMode.WRITE

        print()
        writer.write(traj)
        assert writer.mode == FileWritingMode.APPEND

        captured = capsys.readouterr()
        assert captured.out == """
2 10 10 10 90 90 90

h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
2 11 10 10 90 90 90

h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
"""

        writer = TrajectoryWriter()
        print()
        writer.write(frame1)

        captured = capsys.readouterr()
        assert captured.out == """
2 10 10 10 90 90 90

h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
"""

        frame1 = AtomicSystem(
            atoms=atoms,
            vel=coordinates1,
            cell=Cell(10,
            10,
            10)
        )
        frame2 = AtomicSystem(
            atoms=atoms,
            vel=coordinates2,
            cell=Cell(11,
            10,
            10)
        )

        traj = Trajectory([frame1, frame2])
        writer = TrajectoryWriter()

        print()
        writer.write(traj, traj_type="vel")

        captured = capsys.readouterr()
        assert captured.out == """
2 10 10 10 90 90 90

h 0.000000000000e+00 0.000000000000e+00 0.000000000000e+00
o 0.000000000000e+00 0.000000000000e+00 1.000000000000e+00
2 11 10 10 90 90 90

h 0.000000000000e+00 0.000000000000e+00 0.000000000000e+00
o 0.000000000000e+00 0.000000000000e+00 1.000000000000e+00
"""

        frame1 = AtomicSystem(
            atoms=atoms,
            forces=coordinates1,
            cell=Cell(10,
            10,
            10)
        )
        frame2 = AtomicSystem(
            atoms=atoms,
            forces=coordinates2,
            cell=Cell(11,
            10,
            10)
        )

        traj = Trajectory([frame1, frame2])
        writer = TrajectoryWriter()

        print()
        writer.write(traj, traj_type="force")

        captured = capsys.readouterr()
        assert captured.out == """
2 10 10 10 90 90 90
sum of forces: 0.000000e+00 0.000000e+00 1.000000e+00
h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
2 11 10 10 90 90 90
sum of forces: 0.000000e+00 0.000000e+00 1.000000e+00
h     0.0000000000     0.0000000000     0.0000000000
o     0.0000000000     0.0000000000     1.0000000000
"""

        charges1 = np.array([1, 2])
        charges2 = np.array([3, 4])

        frame1 = AtomicSystem(
            atoms=atoms,
            charges=charges1,
            cell=Cell(10,
            10,
            10)
        )
        frame2 = AtomicSystem(
            atoms=atoms,
            charges=charges2,
            cell=Cell(11,
            10,
            10)
        )

        traj = Trajectory([frame1, frame2])
        writer = TrajectoryWriter()

        print()
        writer.write(traj, traj_type="charge")

        captured = capsys.readouterr()
        assert captured.out == """
2 10 10 10 90 90 90

h 1
o 2
2 11 10 10 90 90 90

h 3
o 4
"""
