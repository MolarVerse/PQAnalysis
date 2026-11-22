import sys

import numpy as np
import pytest

from PQAnalysis.core import Atom, Cell
from PQAnalysis.io import AseTrajReader, TrajectoryReader
from PQAnalysis.io.traj_file.exceptions import AseTrajReaderError
from PQAnalysis.topology import Selection, Topology
from PQAnalysis.traj import Trajectory, TrajectoryFormat

from . import pytestmark  # pylint: disable=unused-import

ase = pytest.importorskip("ase")
# pylint: disable=wrong-import-position,wrong-import-order
from ase.calculators.singlepoint import SinglePointCalculator  # noqa: E402
from ase.io import Trajectory as AseTrajectory  # noqa: E402
import ase.units  # noqa: E402

N_FRAMES = 2
NUMBERS = [8, 1, 1]
# rows are the box vectors, like ASE stores them
BOX = np.array([[10.0, 0.0, 0.0], [2.0, 12.0, 0.0], [0.0, 0.0, 14.0]])
# 1 eV in kcal mol-1, independent of the ASE constants
EV_IN_KCAL_PER_MOL = 23.060548



def positions(frame, n_atoms=len(NUMBERS)):
    return np.arange(3 * n_atoms).reshape(n_atoms, 3) * 0.5 + frame



def velocities_fs(frame, n_atoms=len(NUMBERS)):
    """Velocities in Angstrom fs-1."""
    return np.full((n_atoms, 3), 0.25 * (frame + 1))



def forces_ev(frame, n_atoms=len(NUMBERS)):
    """Forces in eV Angstrom-1."""
    return -positions(frame, n_atoms)



def charges(frame, n_atoms=len(NUMBERS)):
    return np.linspace(-1.0, 1.0, n_atoms) * (frame + 1)



def energy_ev(frame):
    return frame - 5.0



def write_traj(
    filename,
    n_frames=N_FRAMES,
    numbers=NUMBERS,
    pbc=True,
    velocities=True,
    calculator=True,
    initial_charges=False,
    first_step=0,
):
    """Writes a small ASE trajectory shaped like the PQ output."""
    n_atoms = len(numbers)

    with AseTrajectory(filename, "w") as trajectory:
        for frame in range(n_frames):
            atoms = ase.Atoms(
                numbers=numbers,
                positions=positions(frame, n_atoms),
                cell=BOX * (1 + 0.1 * frame),
                pbc=pbc,
            )
            if velocities:
                atoms.set_velocities(
                    velocities_fs(frame, n_atoms) / ase.units.fs
                )
            if initial_charges:
                atoms.set_initial_charges(charges(frame, n_atoms))
            if calculator:
                atoms.calc = SinglePointCalculator(
                    atoms,
                    energy=energy_ev(frame),
                    forces=forces_ev(frame, n_atoms),
                    charges=charges(frame, n_atoms),
                )
            atoms.info["step"] = first_step + 2 * frame
            trajectory.write(atoms)

    return filename



@pytest.fixture
def traj_file(tmp_path):
    return write_traj(str(tmp_path / "traj.traj"))



def test_read_all_properties(traj_file):
    traj = TrajectoryReader(traj_file, traj_format="ase").read()

    assert isinstance(traj, Trajectory)
    assert len(traj) == N_FRAMES

    for i, frame in enumerate(traj):
        assert np.allclose(frame.pos, positions(i))
        assert np.allclose(frame.vel, velocities_fs(i) * 1e15, rtol=1e-12)
        assert np.allclose(
            frame.forces, forces_ev(i) * EV_IN_KCAL_PER_MOL, rtol=1e-6
        )
        assert np.allclose(frame.charges, charges(i))
        assert frame.energy == pytest.approx(
            energy_ev(i) * EV_IN_KCAL_PER_MOL, rel=1e-6
        )

        box = BOX * (1 + 0.1 * i)
        assert np.allclose(frame.cell.box_lengths, np.linalg.norm(box, axis=1))
        assert frame.cell.box_angles[0] == pytest.approx(90.0)
        assert frame.cell.box_angles[1] == pytest.approx(90.0)
        gamma = np.rad2deg(np.arccos(2.0 / np.sqrt(148.0)))
        assert frame.cell.box_angles[2] == pytest.approx(gamma)
        assert np.allclose(frame.cell.volume, abs(np.linalg.det(box)))



def test_topology_from_numbers(tmp_path):
    filename = write_traj(str(tmp_path / "traj.traj"), numbers=[8, 1, 6, 0])

    topology = TrajectoryReader(filename).read()[0].topology
    atoms = topology.atoms

    assert [atom.name for atom in atoms] == ["O", "H", "C", "X"]
    assert [atom.atomic_number for atom in atoms] == [8, 1, 6, 999]
    assert atoms[0] == Atom("O")
    assert np.all(Selection("elem(H)").select(topology) == np.array([1]))



def test_given_topology(traj_file):
    topology = Topology(atoms=[Atom(name) for name in ["OW", "H", "H"]])
    reader = TrajectoryReader(traj_file, topology=topology)

    frames = list(reader.frame_generator())

    assert all(frame.topology is topology for frame in frames)



def test_without_velocities_and_calculator(tmp_path):
    filename = write_traj(
        str(tmp_path / "traj.traj"),
        velocities=False,
        calculator=False,
        initial_charges=True,
    )

    frame = TrajectoryReader(filename).read()[1]

    assert frame.vel.size == 0
    assert frame.forces.size == 0
    assert frame.energy is None
    assert np.allclose(frame.charges, charges(1))



def test_vacuum_cell(tmp_path):
    filename = write_traj(str(tmp_path / "traj.traj"), pbc=False)

    assert TrajectoryReader(filename).read()[0].cell == Cell()



def test_partly_periodic_cell(tmp_path):
    filename = write_traj(
        str(tmp_path / "traj.traj"), pbc=[True, True, False]
    )

    with pytest.raises(AseTrajReaderError) as exception:
        TrajectoryReader(filename).read()

    assert "not periodic in all dimensions" in str(exception.value)



def test_cells(traj_file):
    cells = TrajectoryReader(traj_file).cells

    assert len(cells) == N_FRAMES
    assert np.allclose(cells[1].box_lengths, np.linalg.norm(BOX * 1.1, axis=1))



def test_infer_format(traj_file):
    assert TrajectoryFormat.infer_format_from_extension(traj_file) is \
        TrajectoryFormat.ASE
    assert TrajectoryReader(traj_file).traj_format is TrajectoryFormat.ASE
    assert TrajectoryFormat("ase") is TrajectoryFormat.ASE



def test_multiple_files_and_generators(tmp_path):
    first = write_traj(str(tmp_path / "first.traj"), n_frames=3)
    second = write_traj(str(tmp_path / "second.traj"))

    reader = TrajectoryReader([first, second])

    assert reader.calculate_number_of_frames_per_file() == [3, N_FRAMES]
    assert len(reader.read()) == 3 + N_FRAMES

    frames = list(
        reader.frame_generator(trajectory_start=2, trajectory_stop=4)
    )
    assert len(frames) == 2
    assert np.allclose(frames[0].pos, positions(2))
    assert np.allclose(frames[1].pos, positions(0))

    frames = list(reader.frame_generator(trajectory_start=4))
    assert len(frames) == 1
    assert np.allclose(frames[0].pos, positions(1))



def test_reader_directly(tmp_path):
    filename = write_traj(str(tmp_path / "traj.traj"), first_step=10)
    reader = AseTrajReader(filename)

    assert reader.number_of_frames == N_FRAMES
    assert reader.steps == [10, 12]
    frames = list(reader.frame_generator(start=1))
    assert len(frames) == 1
    assert np.allclose(frames[0].pos, positions(1))



def test_invalid_file(tmp_path):
    filename = tmp_path / "traj.traj"
    filename.write_text("no ase trajectory")

    with pytest.raises(AseTrajReaderError) as exception:
        TrajectoryReader(str(filename)).read()

    assert "Could not open" in str(exception.value)



def test_missing_ase(traj_file, monkeypatch):
    monkeypatch.setitem(sys.modules, "ase", None)

    with pytest.raises(AseTrajReaderError) as exception:
        TrajectoryReader(traj_file).read()

    assert "PQAnalysis[ase]" in str(exception.value)
