"""
Tests of the VACF of the centers of mass of molecules (target_com).

The velocity of every molecule rotates in the xy plane with a common frequency
(random amplitude and phase per molecule), and every atom adds an internal
motion with zero total momentum, so the velocity of a center (the mass
weighted mean of its atoms) is the rotating vector exactly. The normalized
VACF of the centers is then exactly cos(omega * lag), independent of the
time origins, while the VACF of the atoms is not. The results are computed
from the in-memory trajectory and from vel files through the batch kernel, the
streaming kernel (compiled and fallback) and the FFT method.
"""

import sys

import numpy as np
import pytest

from PQAnalysis.analysis import VACF
from PQAnalysis.analysis.vacf import _vacf_kernel_py
from PQAnalysis.analysis.vacf.exceptions import VACFError
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell, Element, Residue
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.topology import Topology
from PQAnalysis.traj import Trajectory

from .. import pytestmark  # pylint: disable=unused-import

# pylint: disable=protected-access

try:
    from PQAnalysis.analysis.vacf import _vacf_kernel
except ModuleNotFoundError:  # pragma: no cover - build-dependent
    _vacf_kernel = None

vacf_module = sys.modules[VACF.__module__]

MOLECULES = [
    ("WAT", ("O", "H", "H"), 6),
    ("MET", ("C", "H", "H", "H", "H"), 3),
    ("ION", ("Na",), 2),
]
N_MOLECULES = sum(count for *_, count in MOLECULES)
N_FRAMES = 40
WINDOW = 10
GAP = 2
OMEGA = 0.3
TIME_STEP = 0.5


def build_topology():
    residues, atoms, residue_ids = [], [], []

    for index, (name, symbols, count) in enumerate(MOLECULES):
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


def rotating_velocities(topology, seed=3):
    """
    The velocities of all frames, shape (n_frames, n_atoms, 3), rounded to
    the digits of the vel file, with the rotating molecule velocities.
    """

    rng = np.random.default_rng(seed)
    masses = np.array([atom.mass for atom in topology.atoms])
    molecule = topology.residue_numbers
    amplitude = rng.uniform(0.5, 2.0, N_MOLECULES)
    phase = rng.uniform(0.0, 2.0 * np.pi, N_MOLECULES)

    frames = []
    for frame in range(N_FRAMES):
        angle = OMEGA * frame + phase
        center = np.stack(
            [amplitude * np.cos(angle), amplitude * np.sin(angle),
             np.zeros(N_MOLECULES)],
            axis=1,
        )
        internal = rng.normal(size=(len(masses), 3))
        for molecule_id in range(N_MOLECULES):
            atoms = molecule == molecule_id
            internal[atoms] -= (
                (internal[atoms] * masses[atoms, None]).sum(0) /
                masses[atoms].sum()
            )
        frames.append(np.round(center[molecule] + internal, 7))

    return np.array(frames)


def expected_vacf():
    return np.cos(OMEGA * np.arange(WINDOW + 1))


def in_memory(velocities, topology):
    return Trajectory([
        AtomicSystem(topology=topology, vel=vel, cell=Cell())
        for vel in velocities
    ])


def write_vel(path, velocities, topology):
    names = [atom.name for atom in topology.atoms]
    lines = []
    for frame in velocities:
        lines.append(f"{len(names)}\n\n")
        lines += [
            f"{name} {x:.7f} {y:.7f} {z:.7f}\n"
            for name, (x, y, z) in zip(names, frame)
        ]
    path.write_text("".join(lines), encoding="utf-8")

    return str(path)


def run_vacf(traj, target="*", **kwargs):
    analysis = VACF(
        traj,
        time_step=TIME_STEP,
        window_size=WINDOW,
        gap=GAP,
        target_species=target,
        target_com=True,
        **kwargs,
    )
    time, vacf = analysis.run()

    return analysis, time, vacf


@pytest.fixture(name="data")
def data_fixture():
    topology = build_topology()

    return topology, rotating_velocities(topology)


class TestInMemory:

    def test_matches_the_analytic_vacf(self, data):
        topology, velocities = data

        analysis, time, vacf = run_vacf(in_memory(velocities, topology))

        assert analysis._raw_reader is None
        assert analysis.n_targets == N_MOLECULES
        assert np.allclose(time, np.arange(WINDOW + 1) * TIME_STEP)
        assert np.allclose(vacf, expected_vacf(), atol=1e-6)

    def test_the_vacf_of_the_atoms_is_not_the_one_of_the_centers(self, data):
        topology, velocities = data

        atoms = VACF(
            in_memory(velocities, topology), time_step=TIME_STEP,
            window_size=WINDOW, gap=GAP,
        ).run()[1]

        assert not np.allclose(atoms, expected_vacf(), atol=1e-2)

    def test_the_fft_method(self, data):
        topology, velocities = data

        analysis, _, vacf = run_vacf(
            in_memory(velocities, topology), method="fft"
        )

        assert analysis.method == "fft"
        # the same estimator for atoms of known VACF: compare with the
        # fft VACF of the centers' own velocities
        masses = np.array([atom.mass for atom in topology.atoms])
        molecule = topology.residue_numbers
        centers = np.array([
            [(frame[molecule == m] * masses[molecule == m, None]).sum(0)
             / masses[molecule == m].sum() for m in range(N_MOLECULES)]
            for frame in velocities
        ])
        reference = VACF(
            Trajectory([
                AtomicSystem(
                    atoms=[Atom("X", use_guess_element=False)] * N_MOLECULES,
                    vel=vel,
                ) for vel in centers
            ]),
            time_step=TIME_STEP, window_size=WINDOW, method="fft",
        ).run()[1]

        assert np.allclose(vacf, reference, atol=1e-10)

    def test_a_selection_defines_the_centers(self, data):
        # the hydrogens of the waters and methanes (their centers do not
        # rotate like the molecule, so only the number of centers is checked)
        topology, velocities = data

        analysis, _, vacf = run_vacf(
            in_memory(velocities, topology), "H"
        )

        assert analysis.n_targets == 9
        assert vacf[0] == pytest.approx(1.0)

    def test_the_default_is_unchanged(self, data):
        topology, velocities = data

        plain = VACF(
            in_memory(velocities, topology), time_step=TIME_STEP,
            window_size=WINDOW, gap=GAP,
        ).run()
        explicit = VACF(
            in_memory(velocities, topology), time_step=TIME_STEP,
            window_size=WINDOW, gap=GAP, target_com=False,
        ).run()

        assert np.array_equal(plain[1], explicit[1])


FILE_PATHS = [
    pytest.param(("python", "stream"), id="fallback-streaming"),
    pytest.param(
        ("compiled", "batch"), id="compiled-batch",
        marks=pytest.mark.skipif(
            _vacf_kernel is None, reason="Cython _vacf_kernel not built"
        ),
    ),
    pytest.param(
        ("compiled", "stream"), id="compiled-streaming",
        marks=pytest.mark.skipif(
            _vacf_kernel is None, reason="Cython _vacf_kernel not built"
        ),
    ),
]


class TestFiles:

    @pytest.fixture(autouse=True)
    def configure(self, request, monkeypatch):
        kernel, path = request.param
        if kernel == "python":
            monkeypatch.setattr(
                vacf_module, "accumulate_frame",
                _vacf_kernel_py.accumulate_frame,
            )
            monkeypatch.setattr(
                vacf_module, "weight_frame", _vacf_kernel_py.weight_frame
            )
        if path == "stream":
            monkeypatch.setattr(VACF, "_direct_batch_max_bytes", 0)

        self.path = path

    @pytest.mark.parametrize("configure", FILE_PATHS, indirect=True)
    def test_matches_the_analytic_vacf(self, configure, data, tmp_path):
        topology, velocities = data
        filename = write_vel(tmp_path / "traj.vel", velocities, topology)

        analysis, _, vacf = run_vacf(
            TrajectoryReader(filename, topology=topology)
        )

        assert analysis._raw_reader is not None
        assert analysis.n_targets == N_MOLECULES
        assert np.allclose(vacf, expected_vacf(), atol=1e-6)

    @pytest.mark.parametrize("configure", FILE_PATHS, indirect=True)
    def test_equals_the_in_memory_trajectory(self, configure, data, tmp_path):
        topology, velocities = data
        filename = write_vel(tmp_path / "traj.vel", velocities, topology)

        from_file = run_vacf(TrajectoryReader(filename, topology=topology))
        from_memory = run_vacf(in_memory(velocities, topology))

        assert np.allclose(from_file[2], from_memory[2], atol=1e-10)

    @pytest.mark.parametrize("configure", FILE_PATHS, indirect=True)
    def test_a_selection_defines_the_centers(self, configure, data, tmp_path):
        topology, velocities = data
        filename = write_vel(tmp_path / "traj.vel", velocities, topology)

        from_file = run_vacf(
            TrajectoryReader(filename, topology=topology), "H"
        )
        from_memory = run_vacf(in_memory(velocities, topology), "H")

        assert from_file[0].n_targets == 9
        assert np.allclose(from_file[2], from_memory[2], atol=1e-10)


class TestErrors:

    @staticmethod
    def traj(topology=None):
        topology = topology or build_topology()
        velocities = rotating_velocities(topology)

        return in_memory(velocities, topology)

    @pytest.mark.parametrize(
        "extra",
        [
            {"charges": np.ones(6 * 3 + 3 * 5 + 2)},
        ],
    )
    def test_the_charge_flux_mode_is_not_supported(self, extra):
        with pytest.raises(VACFError, match="cannot be combined"):
            VACF(
                self.traj(), time_step=TIME_STEP, window_size=WINDOW,
                gap=GAP, target_com=True, **extra,
            )

    def test_a_charge_trajectory_is_not_supported(self):
        topology = build_topology()
        charges = Trajectory([
            AtomicSystem(topology=topology, charges=np.ones(35))
            for _ in range(N_FRAMES)
        ])

        with pytest.raises(VACFError, match="cannot be combined"):
            VACF(
                self.traj(topology), time_step=TIME_STEP,
                window_size=WINDOW, gap=GAP, target_com=True,
                charge_traj=charges,
            )

    def test_a_topology_without_residues_raises(self):
        atoms = [Atom(symbol) for symbol in "OHH" * 4]
        frames = [
            AtomicSystem(atoms=atoms, vel=np.ones((12, 3)) * (f + 1))
            for f in range(N_FRAMES)
        ]

        with pytest.raises(VACFError, match="topology has no residues"):
            VACF(
                Trajectory(frames), time_step=TIME_STEP, window_size=WINDOW,
                gap=GAP, target_com=True,
            )

    def test_atoms_without_a_mass_raise(self):
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
            AtomicSystem(topology=topology, vel=np.ones((4, 3)) * (f + 1))
            for f in range(N_FRAMES)
        ]

        with pytest.raises(VACFError, match="no mass"):
            VACF(
                Trajectory(frames), time_step=TIME_STEP, window_size=WINDOW,
                gap=GAP, target_com=True,
            )
