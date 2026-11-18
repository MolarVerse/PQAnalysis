import sys

import numpy as np
import pytest

from PQAnalysis.core import Atom
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.io.traj_file.exceptions import H5MDReaderError
from PQAnalysis.io.traj_file.h5md_reader import H5MDReader
from PQAnalysis.topology import Selection, Topology
from PQAnalysis.traj import Trajectory, TrajectoryFormat

from . import pytestmark  # pylint: disable=unused-import

h5py = pytest.importorskip("h5py")

N_FRAMES = 3
NAMES = ["O", "H", "H", "X"]
SPECIES = [8, 1, 1, 0]
# rows are the box vectors, like PQ writes them
BOX = np.array([[10.0, 0.0, 0.0], [2.0, 12.0, 0.0], [0.0, 0.0, 14.0]])

UNITS = {
    "position": "Angstrom",
    "velocity": "Angstrom fs-1",
    "force": "kcal mol-1 Angstrom-1",
    "charge": "e",
}



def positions(frame, n_atoms=len(NAMES)):
    return np.arange(3 * n_atoms).reshape(n_atoms, 3) * 0.5 + frame



def velocities(frame, n_atoms=len(NAMES)):
    return np.full((n_atoms, 3), 0.25 * (frame + 1))



def forces(frame, n_atoms=len(NAMES)):
    return -positions(frame, n_atoms)



def charges(frame, n_atoms=len(NAMES)):
    return np.linspace(-1.0, 1.0, n_atoms) * (frame + 1)



def write_h5md(
    filename,
    n_frames=N_FRAMES,
    names=NAMES,
    species=SPECIES,
    elements=("position", "velocity", "force", "charge"),
    units=None,
    box="matrix",
    energy=True,
    group="all",
    first_step=0,
):
    """Writes a small H5MD file shaped like the PQ output."""
    units = {**UNITS, **(units or {})}
    n_atoms = len(NAMES)
    if names is not None or species is not None:
        n_atoms = len(names if names is not None else species)
    steps = np.arange(first_step, first_step + n_frames, dtype=np.int64)
    data = {
        "position": positions,
        "velocity": velocities,
        "force": forces,
        "charge": charges,
    }

    def element(parent, name, values, unit, step=steps):
        element_group = parent.create_group(name)
        element_group.create_dataset("step", data=step)
        element_group.create_dataset("time", data=step * 0.5)
        dataset = element_group.create_dataset("value", data=values)
        if unit is not None:
            dataset.attrs["unit"] = unit
        return element_group

    with h5py.File(filename, "w") as file:
        h5md = file.create_group("h5md")
        h5md.attrs["version"] = np.array([1, 1], dtype=np.int32)

        particles = file.create_group(f"particles/{group}")
        if species is not None:
            particles.create_dataset(
                "species", data=np.array(species, dtype=np.int32)
            )

        if box is not None:
            box_group = particles.create_group("box")
            box_group.attrs["dimension"] = 3
            box_group.attrs["boundary"] = ["periodic"] * 3

        if box == "matrix":
            element(
                box_group,
                "edges",
                np.array([BOX * (1 + 0.1 * i) for i in range(n_frames)]),
                units.get("box", "Angstrom"),
            )
        elif box == "constant":
            edges = box_group.create_dataset("edges", data=np.diag(BOX))
            edges.attrs["unit"] = units.get("box", "Angstrom")

        for name in elements:
            values = np.array(
                [data[name](i, n_atoms) for i in range(n_frames)],
                dtype=np.float32,
            )
            element(particles, name, values, units[name])

        if energy:
            element(
                file.create_group("observables"),
                "potential_energy",
                np.arange(n_frames, dtype=np.float64) - 5.0,
                units.get("energy", "kcal mol-1"),
            )

        if names is not None:
            file.create_dataset(
                "parameters/atom_names",
                data=names,
                dtype=h5py.string_dtype(),
            )

    return filename



@pytest.fixture
def h5md_file(tmp_path):
    return write_h5md(str(tmp_path / "traj.h5md"))



def test_read_all_properties(h5md_file):
    traj = TrajectoryReader(h5md_file, traj_format="h5md").read()

    assert isinstance(traj, Trajectory)
    assert len(traj) == N_FRAMES

    for i, frame in enumerate(traj):
        assert np.allclose(frame.pos, positions(i))
        assert np.allclose(frame.vel, velocities(i) * 1e15)
        assert np.allclose(frame.forces, forces(i))
        assert np.allclose(frame.charges, charges(i), atol=1e-6)
        assert frame.energy == pytest.approx(i - 5.0)

        box = BOX * (1 + 0.1 * i)
        assert np.allclose(frame.cell.box_lengths, np.linalg.norm(box, axis=1))
        assert frame.cell.box_angles[0] == pytest.approx(90.0)
        assert frame.cell.box_angles[1] == pytest.approx(90.0)
        gamma = np.rad2deg(np.arccos(2.0 / np.sqrt(148.0)))
        assert frame.cell.box_angles[2] == pytest.approx(gamma)
        assert np.allclose(frame.cell.volume, abs(np.linalg.det(box)))



def test_topology_from_names_and_species(tmp_path):
    filename = write_h5md(
        str(tmp_path / "traj.h5md"),
        names=["OW", "Co", "h", "X", "dum"],
        species=[8, 6, 1, 0, 0],
    )

    topology = TrajectoryReader(filename).read()[0].topology
    atoms = topology.atoms

    assert [atom.name for atom in atoms] == ["OW", "Co", "h", "X", "dum"]
    assert np.all(Selection("elem(O)").select(topology) == np.array([0]))
    assert atoms[0] == Atom("OW", "o")
    assert atoms[1] == Atom("Co", "c")
    assert atoms[2] == Atom("h")
    assert atoms[3] == Atom("X")
    assert atoms[4] == Atom("dum")



def test_topology_like_xyz_without_species(tmp_path):
    filename = write_h5md(
        str(tmp_path / "traj.h5md"),
        names=["O", "OW"],
        species=None,
    )

    atoms = TrajectoryReader(filename).read()[0].topology.atoms

    assert atoms[0] == Atom("O")
    assert atoms[1] == Atom("OW", use_guess_element=False)



def test_topology_from_species_without_names(tmp_path):
    filename = write_h5md(str(tmp_path / "traj.h5md"), names=None)

    atoms = TrajectoryReader(filename).read()[0].topology.atoms

    assert [atom.name for atom in atoms] == ["o", "h", "h", "X"]
    assert [atom.atomic_number for atom in atoms] == [8, 1, 1, 999]



def test_given_topology_and_selection(h5md_file):
    topology = Topology(atoms=[Atom(name) for name in ["O", "H", "H", "X"]])
    reader = TrajectoryReader(h5md_file, topology=topology)

    frames = list(reader.frame_generator())

    assert all(frame.topology is topology for frame in frames)
    assert np.all(
        Selection("elem(H)").select(frames[0].topology) == np.array([1, 2])
    )



@pytest.mark.parametrize(
    "unit, factor",
    [
        ("Angstrom fs-1", 1e15),
        ("Angstrom ps-1", 1e12),
        ("Angstrom s-1", 1.0),
        ("nm ps-1", 1e13),
        (b"Angstrom  fs-1", 1e15),
        (np.array([b"Angstrom fs-1"]), 1e15),
    ],
)
def test_velocity_units(tmp_path, unit, factor):
    filename = write_h5md(
        str(tmp_path / "traj.h5md"), units={"velocity": unit}
    )

    frame = TrajectoryReader(filename).read()[1]

    assert np.allclose(frame.vel, velocities(1) * factor)



def test_length_units_nm(tmp_path):
    filename = write_h5md(
        str(tmp_path / "traj.h5md"),
        units={
            "position": "nm",
            "box": "nm",
        },
    )

    frame = TrajectoryReader(filename).read()[0]

    assert np.allclose(frame.pos, positions(0) * 10.0)
    assert np.allclose(
        frame.cell.box_lengths, np.linalg.norm(BOX, axis=1) * 10.0
    )



@pytest.mark.parametrize(
    "units",
    [
        {
            "velocity": "m s-1"
        },
        {
            "velocity": None
        },
        {
            "velocity": np.array([b"Angstrom", b"fs-1"])
        },
        {
            "position": "bohr"
        },
        {
            "force": "eV Angstrom-1"
        },
        {
            "charge": "C"
        },
        {
            "energy": "eV"
        },
        {
            "box": "pm"
        },
    ],
)
def test_unsupported_unit(tmp_path, units):
    filename = write_h5md(str(tmp_path / "traj.h5md"), units=units)

    with pytest.raises(H5MDReaderError) as exception:
        TrajectoryReader(filename).read()

    assert "is not supported" in str(exception.value)



def test_missing_datasets(tmp_path):
    filename = write_h5md(
        str(tmp_path / "traj.h5md"),
        elements=("velocity", ),
        box=None,
        energy=False,
    )

    reader = TrajectoryReader(filename)
    frame = reader.read()[2]

    assert np.allclose(frame.vel, velocities(2) * 1e15)
    assert frame.pos.shape == (0, 3)
    assert frame.forces.shape == (0, 3)
    assert frame.energy is None
    assert frame.cell.is_vacuum
    assert all(cell.is_vacuum for cell in reader.cells)



def test_constant_box(tmp_path):
    filename = write_h5md(str(tmp_path / "traj.h5md"), box="constant")

    reader = TrajectoryReader(filename)

    assert len(reader.cells) == N_FRAMES
    for cell in reader.cells:
        assert np.allclose(cell.box_lengths, [10.0, 12.0, 14.0])
        assert np.allclose(cell.box_angles, 90.0)



def test_cells(h5md_file):
    cells = TrajectoryReader(h5md_file).cells

    assert len(cells) == N_FRAMES
    assert np.allclose(cells[2].box_lengths, np.linalg.norm(BOX * 1.2, axis=1))



@pytest.mark.parametrize("extension", [".h5md", ".h5"])
def test_infer_format(tmp_path, extension):
    filename = write_h5md(str(tmp_path / f"traj{extension}"))

    assert TrajectoryFormat.infer_format_from_extension(filename) is \
        TrajectoryFormat.H5MD
    assert TrajectoryReader(filename).traj_format is TrajectoryFormat.H5MD
    assert TrajectoryFormat("h5md") is TrajectoryFormat.H5MD



def test_multiple_files_and_generators(tmp_path):
    first = write_h5md(str(tmp_path / "first.h5md"))
    second = write_h5md(str(tmp_path / "second.h5md"), n_frames=2)

    reader = TrajectoryReader([first, second])

    assert reader.calculate_number_of_frames_per_file() == [N_FRAMES, 2]
    assert len(reader.read()) == N_FRAMES + 2

    frames = list(
        reader.frame_generator(trajectory_start=2, trajectory_stop=4)
    )
    assert len(frames) == 2
    assert np.allclose(frames[0].pos, positions(2))
    assert np.allclose(frames[1].pos, positions(0))

    frames = list(reader.frame_generator(trajectory_start=4))
    assert len(frames) == 1
    assert np.allclose(frames[0].pos, positions(1))

    windows = list(reader.window_generator(window_size=2, window_gap=3))
    assert len(windows) == 2
    assert np.allclose(windows[1][0].pos, positions(0))
    assert np.allclose(windows[1][1].pos, positions(1))



def test_h5md_reader_directly(h5md_file):
    reader = H5MDReader(h5md_file)

    assert reader.number_of_frames == N_FRAMES
    frames = list(reader.frame_generator(start=1))
    assert len(frames) == 2
    assert np.allclose(frames[0].pos, positions(1))



def test_other_particle_group(tmp_path):
    filename = write_h5md(str(tmp_path / "traj.h5md"), group="water")

    assert len(TrajectoryReader(filename).read()) == N_FRAMES



def test_several_particle_groups(tmp_path):
    filename = write_h5md(str(tmp_path / "traj.h5md"), group="water")
    with h5py.File(filename, "a") as file:
        file.create_group("particles/ions")

    with pytest.raises(H5MDReaderError) as exception:
        TrajectoryReader(filename).read()

    assert "single particle group" in str(exception.value)



def test_different_steps(tmp_path):
    filename = write_h5md(str(tmp_path / "traj.h5md"))
    with h5py.File(filename, "a") as file:
        file["particles/all/velocity/step"][0] = 7

    with pytest.raises(H5MDReaderError) as exception:
        TrajectoryReader(filename).read()

    assert "other steps" in str(exception.value)



def test_no_particle_data(tmp_path):
    filename = write_h5md(str(tmp_path / "traj.h5md"), elements=())

    with pytest.raises(H5MDReaderError) as exception:
        TrajectoryReader(filename).read()

    assert "no position, velocity, force or charge" in str(exception.value)



def test_no_names_and_species(tmp_path):
    filename = write_h5md(
        str(tmp_path / "traj.h5md"), names=None, species=None
    )

    with pytest.raises(H5MDReaderError) as exception:
        TrajectoryReader(filename).read()

    assert "neither" in str(exception.value)



def test_time_dependent_species(tmp_path):
    filename = write_h5md(str(tmp_path / "traj.h5md"), species=None)
    with h5py.File(filename, "a") as file:
        file.create_group("particles/all/species")

    with pytest.raises(H5MDReaderError) as exception:
        TrajectoryReader(filename).read()

    assert "Time-dependent species" in str(exception.value)



def test_non_periodic_box(tmp_path):
    filename = write_h5md(str(tmp_path / "traj.h5md"))
    with h5py.File(filename, "a") as file:
        file["particles/all/box"].attrs["boundary"] = [
            "periodic",
            "none",
            "periodic",
        ]

    with pytest.raises(H5MDReaderError) as exception:
        TrajectoryReader(filename).read()

    assert "not periodic" in str(exception.value)



def test_invalid_file(tmp_path):
    filename = tmp_path / "traj.h5md"
    filename.write_text("no hdf5")

    with pytest.raises(H5MDReaderError) as exception:
        TrajectoryReader(str(filename)).read()

    assert "Could not open" in str(exception.value)



def test_missing_h5py(h5md_file, monkeypatch):
    monkeypatch.setitem(sys.modules, "h5py", None)

    with pytest.raises(H5MDReaderError) as exception:
        TrajectoryReader(h5md_file).read()

    assert "PQAnalysis[h5md]" in str(exception.value)
