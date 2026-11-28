import sys

import numpy as np
import pytest

from PQAnalysis.core import Atom
from PQAnalysis.exceptions import PQFileNotFoundError
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.io.traj_file.exceptions import ZarrTrajReaderError
from PQAnalysis.io.traj_file.zarr_traj_reader import ZarrTrajReader
from PQAnalysis.traj import Trajectory, TrajectoryFormat

from . import pytestmark  # pylint: disable=unused-import
from .test_h5mdReader import (
    BOX,
    N_FRAMES,
    NAMES,
    SPECIES,
    UNITS,
    charges,
    forces,
    positions,
    velocities,
    write_h5md,
)

zarr = pytest.importorskip("zarr")



def write_zarr(
    filename,
    n_frames=N_FRAMES,
    names=NAMES,
    species=SPECIES,
    elements=("position", "velocity", "force", "charge"),
    units=None,
    box="matrix",
    energy=True,
    frames_per_chunk=1,
    frames_per_shard=None,
    names_as_array=False,
):
    """Writes a small Zarr store in the H5MD layout of PQ."""
    units = {**UNITS, **(units or {})}
    n_atoms = len(names if names is not None else species)
    steps = np.arange(n_frames, dtype=np.int64)
    data = {
        "position": positions,
        "velocity": velocities,
        "force": forces,
        "charge": charges,
    }

    def element(parent, name, values, unit):
        element_group = parent.create_group(name)
        element_group.create_array("step", data=steps)
        element_group.create_array("time", data=steps * 0.5)
        array = element_group.create_array(
            "value",
            data=values,
            chunks=(frames_per_chunk, ) + values.shape[1:],
            shards=None if frames_per_shard is None else
            (frames_per_shard, ) + values.shape[1:],
        )
        array.attrs["unit"] = unit

    root = zarr.open_group(filename, mode="w")
    root.create_group("h5md").attrs["version"] = [1, 1]

    particles = root.create_group("particles/all")
    if species is not None:
        particles.create_array(
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
        edges = box_group.create_array("edges", data=np.diag(BOX))
        edges.attrs["unit"] = units.get("box", "Angstrom")

    for name in elements:
        values = np.array(
            [data[name](i, n_atoms) for i in range(n_frames)],
            dtype=np.float32,
        )
        element(particles, name, values, units[name])

    if energy:
        element(
            root.create_group("observables"),
            "potential_energy",
            np.arange(n_frames, dtype=np.float64) - 5.0,
            units.get("energy", "kcal mol-1"),
        )

    if names is not None and names_as_array:
        root.create_group("parameters").create_array(
            "atom_names", data=np.array(names, dtype=str)
        )
    elif names is not None:
        root.create_group("parameters").attrs["atom_names"] = list(names)

    return filename



@pytest.fixture
def zarr_store(tmp_path):
    return write_zarr(str(tmp_path / "traj.zarr"))



def test_read_all_properties(zarr_store):
    traj = TrajectoryReader(zarr_store, traj_format="zarr").read()

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
        assert np.allclose(frame.cell.volume, abs(np.linalg.det(box)))



@pytest.mark.parametrize("frames_per_chunk", [1, 2, 10])
def test_same_frames_as_h5md(tmp_path, frames_per_chunk):
    store = write_zarr(
        str(tmp_path / "traj.zarr"),
        n_frames=5,
        frames_per_chunk=frames_per_chunk,
    )
    h5md = write_h5md(str(tmp_path / "traj.h5md"), n_frames=5)

    zarr_traj = TrajectoryReader(store).read()
    h5md_traj = TrajectoryReader(h5md).read()

    assert len(zarr_traj) == len(h5md_traj) == 5
    for zarr_frame, h5md_frame in zip(zarr_traj, h5md_traj):
        assert zarr_frame.topology == h5md_frame.topology
        assert zarr_frame.cell == h5md_frame.cell
        assert np.array_equal(zarr_frame.pos, h5md_frame.pos)
        assert np.array_equal(zarr_frame.vel, h5md_frame.vel)
        assert np.array_equal(zarr_frame.forces, h5md_frame.forces)
        assert np.array_equal(zarr_frame.charges, h5md_frame.charges)
        assert zarr_frame.energy == h5md_frame.energy



@pytest.mark.parametrize(
    "frames_per_chunk, frames_per_shard", [
        (3, None),
        (1, 2),
    ]
)
def test_chunked_frames_out_of_order(
    tmp_path, frames_per_chunk, frames_per_shard
):
    store = write_zarr(
        str(tmp_path / "traj.zarr"),
        n_frames=7,
        frames_per_chunk=frames_per_chunk,
        frames_per_shard=frames_per_shard,
    )

    reader = TrajectoryReader(store)
    frames = list(reader.frame_generator(trajectory_start=2))

    assert len(frames) == 5
    for i, frame in enumerate(frames, start=2):
        assert np.allclose(frame.pos, positions(i))

    frames = list(ZarrTrajReader(store).frame_generator(start=4, stop=6))
    assert [frame.energy for frame in frames] == [-1.0, 0.0]



def test_topology_from_names_and_species(tmp_path):
    store = write_zarr(
        str(tmp_path / "traj.zarr"),
        names=["OW", "h", "dum"],
        species=[8, 1, 0],
    )

    atoms = TrajectoryReader(store).read()[0].topology.atoms

    assert atoms[0] == Atom("OW", "o")
    assert atoms[1] == Atom("h")
    assert atoms[2] == Atom("dum")



def test_atom_names_as_array(tmp_path):
    store = write_zarr(
        str(tmp_path / "traj.zarr"),
        names=["OW", "h"],
        species=[8, 1],
        names_as_array=True,
    )

    atoms = TrajectoryReader(store).read()[0].topology.atoms

    assert [atom.name for atom in atoms] == ["OW", "h"]



def test_parameters_without_atom_names(tmp_path):
    store = write_zarr(str(tmp_path / "traj.zarr"), names=None)
    zarr.open_group(store, mode="a").create_group("parameters")

    atoms = TrajectoryReader(store).read()[0].topology.atoms

    assert [atom.name for atom in atoms] == ["o", "h", "h", "X"]



def test_topology_from_species_without_names(tmp_path):
    store = write_zarr(str(tmp_path / "traj.zarr"), names=None)

    atoms = TrajectoryReader(store).read()[0].topology.atoms

    assert [atom.name for atom in atoms] == ["o", "h", "h", "X"]



def test_missing_arrays_and_constant_box(tmp_path):
    store = write_zarr(
        str(tmp_path / "traj.zarr"),
        elements=("velocity", ),
        box="constant",
        energy=False,
    )

    reader = TrajectoryReader(store)
    frame = reader.read()[2]

    assert np.allclose(frame.vel, velocities(2) * 1e15)
    assert frame.pos.shape == (0, 3)
    assert frame.energy is None
    assert len(reader.cells) == N_FRAMES
    assert np.allclose(reader.cells[1].box_lengths, [10.0, 12.0, 14.0])



def test_no_box(tmp_path):
    store = write_zarr(str(tmp_path / "traj.zarr"), box=None)

    assert all(cell.is_vacuum for cell in TrajectoryReader(store).cells)



def test_multiple_stores(tmp_path):
    first = write_zarr(str(tmp_path / "first.zarr"))
    second = write_zarr(str(tmp_path / "second.zarr"), n_frames=2)

    reader = TrajectoryReader([first, second])

    assert reader.calculate_number_of_frames_per_file() == [N_FRAMES, 2]
    frames = list(
        reader.frame_generator(trajectory_start=2, trajectory_stop=4)
    )
    assert np.allclose(frames[0].pos, positions(2))
    assert np.allclose(frames[1].pos, positions(0))



def test_unsupported_unit(tmp_path):
    store = write_zarr(
        str(tmp_path / "traj.zarr"), units={"velocity": "m s-1"}
    )

    with pytest.raises(ZarrTrajReaderError) as exception:
        TrajectoryReader(store).read()

    assert "Zarr dataset /particles/all/velocity/value" in str(exception.value)



def test_time_dependent_species(tmp_path):
    store = write_zarr(str(tmp_path / "traj.zarr"), species=None)
    zarr.open_group(store, mode="a").create_group("particles/all/species")

    with pytest.raises(ZarrTrajReaderError) as exception:
        TrajectoryReader(store).read()

    assert "Time-dependent species" in str(exception.value)



def test_no_particle_data(tmp_path):
    store = write_zarr(str(tmp_path / "traj.zarr"), elements=())

    with pytest.raises(ZarrTrajReaderError) as exception:
        TrajectoryReader(store).read()

    assert "The Zarr file" in str(exception.value)



@pytest.mark.parametrize("suffix", ["", "/"])
def test_infer_format(zarr_store, suffix):
    filename = zarr_store + suffix

    assert TrajectoryFormat.infer_format_from_extension(filename) is \
        TrajectoryFormat.ZARR
    assert TrajectoryReader(filename).traj_format is TrajectoryFormat.ZARR
    assert TrajectoryFormat("zarr") is TrajectoryFormat.ZARR
    assert len(TrajectoryReader(filename).read()) == N_FRAMES



def test_zarr_v2_store(tmp_path):
    store = str(tmp_path / "traj.zarr")
    root = zarr.open_group(store, mode="w", zarr_format=2)
    particles = root.create_group("particles/all")
    position = particles.create_group("position")
    position.create_array("step", data=np.arange(2))
    value = position.create_array(
        "value", data=np.array([positions(0), positions(1)])
    )
    value.attrs["unit"] = "Angstrom"
    particles.create_array("species", data=np.array(SPECIES))

    traj = TrajectoryReader(store).read()

    assert len(traj) == 2
    assert np.allclose(traj[1].pos, positions(1))



def test_directory_without_store(tmp_path):
    directory = tmp_path / "traj.zarr"
    directory.mkdir()

    with pytest.raises(PQFileNotFoundError):
        TrajectoryReader(str(directory))

    with pytest.raises(PQFileNotFoundError):
        TrajectoryReader([str(directory)])



def test_invalid_store(tmp_path):
    directory = tmp_path / "traj.zarr"
    directory.mkdir()
    (directory / "zarr.json").write_text("no json")

    with pytest.raises(ZarrTrajReaderError) as exception:
        TrajectoryReader(str(directory)).read()

    assert "Could not open the Zarr store" in str(exception.value)



def test_missing_zarr(zarr_store, monkeypatch):
    monkeypatch.setitem(sys.modules, "zarr", None)

    with pytest.raises(ZarrTrajReaderError) as exception:
        TrajectoryReader(zarr_store).read()

    assert "PQAnalysis[zarr]" in str(exception.value)
