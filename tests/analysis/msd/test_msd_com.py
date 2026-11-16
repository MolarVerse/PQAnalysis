"""
Tests of the MSD of the centers of mass of molecules (target_com).

Waters, methanes and ions move with a constant velocity and cross the
periodic boundary several times, every atom is wrapped into the box on its
own. The MSD of such a motion is known exactly,

    MSD(lag) = mean over the molecules of (v_i * lag)^2   per axis,

so the results are compared with this analytic value, independent of any
unwrapping or center of mass code. They are computed from the in-memory
trajectory and from files through every code path of the MSD class: the exact
batch kernel, the exact streaming kernel and the bounded streaming kernel
(compiled and fallback).
"""

import sys

import numpy as np
import pytest

from scipy.spatial.transform import Rotation

from PQAnalysis.analysis import MSD
from PQAnalysis.analysis.msd import _msd_kernel_py
from PQAnalysis.analysis.msd.exceptions import MSDError
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.atomic_system.exceptions import AtomicSystemError
from PQAnalysis.core import Atom, Cell, Element, Residue
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.topology import Topology
from PQAnalysis.traj import Trajectory

from .. import pytestmark  # pylint: disable=unused-import

# pylint: disable=protected-access

try:
    from PQAnalysis.analysis.msd import _msd_kernel
except ModuleNotFoundError:  # pragma: no cover - build-dependent
    _msd_kernel = None

msd_module = sys.modules[MSD.__module__]

WATER = np.array([[0.0, 0, 0], [0.9572, 0, 0], [-0.2399, 0.9266, 0]])
METHANE = np.array([
    [0.0, 0, 0],
    [0.63, 0.63, 0.63],
    [-0.63, -0.63, 0.63],
    [-0.63, 0.63, -0.63],
    [0.63, -0.63, -0.63],
])

MOLECULES = [
    ("WAT", ("O", "H", "H"), WATER, 6),
    ("MET", ("C", "H", "H", "H", "H"), METHANE, 3),
    ("ION", ("Na",), np.zeros((1, 3)), 2),
]
N_MOLECULES = sum(count for *_, count in MOLECULES)
N_FRAMES = 30
WINDOW = 10
GAP = 2
CELLS = {
    "orthorhombic": Cell(10.0, 12.0, 9.0),
    "triclinic": Cell(10.0, 11.0, 9.0, 80.0, 95.0, 100.0),
}


def build_topology():
    residues, atoms, residue_ids = [], [], []

    for index, (name, symbols, _, count) in enumerate(MOLECULES):
        residues.append(
            Residue(
                name=name,
                elements=[Element(symbol) for symbol in symbols],
                residue_id=index + 1,
                total_charge=0,
                atom_types=np.arange(len(symbols)),
                partial_charges=np.zeros(len(symbols)),
            )
        )
        for _ in range(count):
            atoms += [Atom(symbol) for symbol in symbols]
            residue_ids += [index + 1] * len(symbols)

    return Topology(
        atoms=atoms,
        residue_ids=np.array(residue_ids),
        reference_residues=residues,
    )


def moving_molecules(cell, seed=11):
    """
    The wrapped frames and the Cartesian velocity (per frame) of every
    molecule. Fractional velocities up to 0.4 box per frame, so the
    molecules cross the boundary every few frames.
    """

    rng = np.random.default_rng(seed)
    matrix = cell.box_matrix
    velocities = rng.uniform(-0.4, 0.4, (N_MOLECULES, 3)) @ matrix.T
    starts = rng.uniform(0.0, 1.0, (N_MOLECULES, 3)) @ matrix.T
    orientations = Rotation.from_quat(rng.normal(size=(N_MOLECULES, 4)))
    topology = build_topology()
    offsets = []
    for _, symbols, geometry, count in MOLECULES:
        masses = np.array([Atom(symbol).mass for symbol in symbols])
        centered = geometry - (masses[:, None] * geometry).sum(0) / masses.sum()
        offsets += [centered] * count
    rotated = [
        orientation.apply(offset)
        for orientation, offset in zip(orientations, offsets)
    ]

    frames = []
    for frame in range(N_FRAMES):
        atoms = np.vstack([
            rotated[i] + starts[i] + frame * velocities[i]
            for i in range(N_MOLECULES)
        ])
        frames.append(
            np.round(cell.image(atoms), 6)
        )

    return frames, velocities, topology


def analytic_msd(velocities, lags):
    """The MSD per axis: the mean over the molecules of (v * lag)^2."""

    return (velocities**2).mean(axis=0)[None, :] * lags[:, None]**2


def in_memory(frames, topology, cell):
    return Trajectory([
        AtomicSystem(topology=topology, pos=positions, cell=cell)
        for positions in frames
    ])


def write_xyz(path, frames, topology, cell):
    names = [atom.name for atom in topology.atoms]
    lengths = " ".join(str(value) for value in cell.box_lengths)
    angles = " ".join(str(value) for value in cell.box_angles)
    lines = []
    for positions in frames:
        lines.append(f"{len(names)} {lengths} {angles}\n\n")
        lines += [
            f"{name} {x:.6f} {y:.6f} {z:.6f}\n"
            for name, (x, y, z) in zip(names, positions)
        ]
    path.write_text("".join(lines), encoding="utf-8")

    return str(path)


def run_msd(traj, target="*", **kwargs):
    analysis = MSD(
        traj, target, window=WINDOW, gap=GAP, target_com=True, **kwargs
    )
    lags, msd_x, msd_y, msd_z, msd_tot = analysis.run()

    return analysis, lags, np.column_stack([msd_x, msd_y, msd_z]), msd_tot


def assert_analytic(lags, msd, msd_tot, velocities):
    expected = analytic_msd(velocities, lags)

    assert np.allclose(msd, expected, rtol=1e-4, atol=1e-6)
    assert np.allclose(msd_tot, expected.sum(axis=1), rtol=1e-4, atol=1e-6)
    assert msd_tot[-1] > 50.0  # many box lengths squared


@pytest.fixture(name="system", params=CELLS.values(), ids=CELLS.keys())
def system_fixture(request):
    cell = request.param
    frames, velocities, topology = moving_molecules(cell)

    return cell, frames, velocities, topology


class TestInMemory:

    def test_matches_the_analytic_msd(self, system):
        cell, frames, velocities, topology = system

        analysis, lags, msd, msd_tot = run_msd(
            in_memory(frames, topology, cell)
        )

        assert analysis._raw_reader is None
        assert len(analysis.target_indices) == N_MOLECULES
        assert_analytic(lags, msd, msd_tot, velocities)

    def test_the_msd_of_atoms_is_not_the_msd_of_molecules(self, system):
        # the oxygen moves like its molecule, but a mean over atoms weights
        # the molecules by their size
        cell, frames, velocities, topology = system

        centers = run_msd(in_memory(frames, topology, cell))[3]
        atoms = MSD(
            in_memory(frames, topology, cell), "*", window=WINDOW, gap=GAP
        ).run()[4]

        assert not np.allclose(centers, atoms, rtol=1e-3)

    def test_a_selection_defines_the_centers(self, system):
        # the hydrogens of the waters and methanes: the ions have none, and
        # the center of the hydrogens moves with its molecule
        cell, frames, velocities, topology = system

        analysis, lags, msd, msd_tot = run_msd(
            in_memory(frames, topology, cell), "H"
        )

        assert len(analysis.target_indices) == 9
        assert_analytic(lags, msd, msd_tot, velocities[:9])

    def test_the_fit_uses_the_molecules(self, system):
        cell, frames, velocities, topology = system

        analysis = MSD(
            in_memory(frames, topology, cell), "*", window=WINDOW, gap=GAP,
            time_step=0.5, fit_window=5, target_com=True,
        )
        analysis.run()

        slope = np.polyfit(
            analysis.lags[-5:] * 0.5, analysis.msd_tot[-5:], 1
        )[0]
        assert analysis.fit_results["total"].slope == pytest.approx(slope)


FILE_PATHS = [
    pytest.param("python", "stream", id="fallback-bounded-streaming"),
    pytest.param(
        "compiled", "batch", id="compiled-batch",
        marks=pytest.mark.skipif(
            _msd_kernel is None, reason="Cython _msd_kernel not built"
        ),
    ),
    pytest.param(
        "compiled", "exact-stream", id="compiled-exact-streaming",
        marks=pytest.mark.skipif(
            _msd_kernel is None, reason="Cython _msd_kernel not built"
        ),
    ),
]


class TestFiles:

    @pytest.fixture(autouse=True)
    def configure(self, request, monkeypatch):
        kernel, path = request.node.callspec.params["kernel_path"]

        if kernel == "python":
            monkeypatch.setattr(
                msd_module, "msd_frame_update",
                _msd_kernel_py.msd_frame_update,
            )
        if path == "exact-stream":
            monkeypatch.setattr(MSD, "_direct_batch_max_bytes", 0)

        self.path = path

    @pytest.mark.parametrize(
        "kernel_path", [(p.values[0], p.values[1]) for p in FILE_PATHS],
        ids=[p.id for p in FILE_PATHS],
    )
    def test_matches_the_analytic_msd(self, kernel_path, system, tmp_path):
        cell, frames, velocities, topology = system
        if cell.alpha != 90.0 and self.path != "stream":
            pytest.skip("the exact kernels are orthorhombic only")
        filename = write_xyz(tmp_path / "traj.xyz", frames, topology, cell)

        analysis, lags, msd, msd_tot = run_msd(
            TrajectoryReader(filename, topology=topology)
        )

        assert analysis._raw_reader is not None
        assert len(analysis.target_indices) == N_MOLECULES
        assert_analytic(lags, msd, msd_tot, velocities)

    @pytest.mark.parametrize(
        "kernel_path", [(p.values[0], p.values[1]) for p in FILE_PATHS],
        ids=[p.id for p in FILE_PATHS],
    )
    def test_equals_the_in_memory_trajectory(
        self, kernel_path, system, tmp_path
    ):
        cell, frames, velocities, topology = system
        if cell.alpha != 90.0 and self.path != "stream":
            pytest.skip("the exact kernels are orthorhombic only")
        filename = write_xyz(tmp_path / "traj.xyz", frames, topology, cell)

        from_file = run_msd(TrajectoryReader(filename, topology=topology))
        from_memory = run_msd(in_memory(frames, topology, cell))

        assert np.allclose(from_file[2], from_memory[2], rtol=1e-10, atol=1e-12)
        assert np.allclose(from_file[3], from_memory[3], rtol=1e-10, atol=1e-12)

    @pytest.mark.parametrize(
        "kernel_path", [(p.values[0], p.values[1]) for p in FILE_PATHS],
        ids=[p.id for p in FILE_PATHS],
    )
    def test_a_selection_defines_the_centers(
        self, kernel_path, system, tmp_path
    ):
        cell, frames, velocities, topology = system
        if cell.alpha != 90.0 and self.path != "stream":
            pytest.skip("the exact kernels are orthorhombic only")
        filename = write_xyz(tmp_path / "traj.xyz", frames, topology, cell)

        analysis, lags, msd, msd_tot = run_msd(
            TrajectoryReader(filename, topology=topology), "H"
        )

        assert len(analysis.target_indices) == 9
        assert_analytic(lags, msd, msd_tot, velocities[:9])


class TestWithoutCenters:

    def test_the_default_is_unchanged(self, system):
        cell, frames, _, topology = system
        traj = in_memory(frames, topology, cell)

        plain = MSD(traj, "O", window=WINDOW, gap=GAP).run()
        explicit = MSD(
            in_memory(frames, topology, cell), "O", window=WINDOW, gap=GAP,
            target_com=False,
        ).run()

        for a, b in zip(plain, explicit):
            assert np.array_equal(a, b)


class TestErrors:

    def test_a_topology_without_residues_raises(self):
        atoms = [Atom(symbol) for symbol in "OHH" * 4]
        cell = CELLS["orthorhombic"]
        frames = [
            AtomicSystem(atoms=atoms, pos=pos[:12], cell=cell)
            for pos in moving_molecules(cell)[0]
        ]

        with pytest.raises(MSDError, match="topology has no residues"):
            MSD(
                Trajectory(frames), "O", window=WINDOW, gap=GAP,
                target_com=True,
            )

    def test_atoms_without_a_mass_raise(self):
        cell = CELLS["orthorhombic"]
        residue = Residue(
            name="MOL", elements=[Element("C"), Element("C")], residue_id=1,
            total_charge=0, atom_types=np.array([0, 0]),
            partial_charges=np.zeros(2),
        )
        topology = Topology(
            atoms=[Atom(n, use_guess_element=False) for n in ("A1", "A2")] * 2,
            residue_ids=np.ones(4, dtype=int),
            reference_residues=[residue],
        )
        frames = [
            AtomicSystem(
                topology=topology,
                pos=np.arange(12, dtype=float).reshape(4, 3) + frame,
                cell=cell,
            ) for frame in range(40)
        ]

        with pytest.raises(MSDError, match="no mass"):
            MSD(
                Trajectory(frames), "*", window=WINDOW, gap=GAP,
                target_com=True,
            )

    def test_a_molecule_larger_than_half_the_box_raises_when_running(self):
        cell = CELLS["orthorhombic"]
        residue = Residue(
            name="LONG", elements=[Element("C")] * 3, residue_id=1,
            total_charge=0, atom_types=np.zeros(3, dtype=int),
            partial_charges=np.zeros(3),
        )
        topology = Topology(
            atoms=[Atom("C")] * 3, residue_ids=np.ones(3, dtype=int),
            reference_residues=[residue],
        )
        frames = [
            AtomicSystem(
                topology=topology,
                pos=np.array([[0.0, 0, 0], [3.2, 0, 0], [6.4, 0, 0]]),
                cell=cell,
            ) for _ in range(40)
        ]
        analysis = MSD(
            Trajectory(frames), "*", window=WINDOW, gap=GAP, target_com=True
        )

        with pytest.raises(AtomicSystemError, match="spans half of the box"):
            analysis.run()

    def test_placeholder_atoms_log_a_warning(self, caplog):
        cell = CELLS["orthorhombic"]
        residue = Residue(
            name="MOL",
            elements=[Element("O"), Element("H"), Element("X")],
            residue_id=1, total_charge=0, atom_types=np.arange(3),
            partial_charges=np.zeros(3),
        )
        topology = Topology(
            atoms=[Atom(name) for name in "OHX"],
            residue_ids=np.ones(3, dtype=int), reference_residues=[residue],
        )
        frames = [
            AtomicSystem(
                topology=topology,
                pos=np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]]) + 0.1 * f,
                cell=cell,
            ) for f in range(40)
        ]

        with caplog.at_level("WARNING"):
            MSD(
                Trajectory(frames), "*", window=WINDOW, gap=GAP,
                target_com=True,
            )

        assert "placeholder mass (x)" in caplog.text.lower()
