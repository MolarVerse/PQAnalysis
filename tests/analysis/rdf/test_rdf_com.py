"""
Tests of the RDF between the centers of mass of molecules (reference_com and
target_com).

A box of waters, methane molecules and ions is generated with seeded random
molecule positions and orientations, and every atom is wrapped into the box on
its own, so that molecules are split by the periodic boundary. The histogram
of the RDF is compared exactly with an independent brute force calculation
(true centers of mass, minimum image distances) for the in-memory trajectory,
the streamed file and the batch read file, with the compiled and the fallback
kernel.
"""

import sys

import numpy as np
import pytest

from scipy.spatial.transform import Rotation

from PQAnalysis.analysis import RDF
from PQAnalysis.analysis.rdf import _rdf_kernel_py
from PQAnalysis.analysis.rdf.exceptions import RDFError
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.atomic_system.exceptions import AtomicSystemError
from PQAnalysis.core import Atom, Cell, Element, Residue
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.topology import Topology
from PQAnalysis.traj import Trajectory

from .. import pytestmark  # pylint: disable=unused-import

# pylint: disable=protected-access

try:
    from PQAnalysis.analysis.rdf import _rdf_kernel
except ModuleNotFoundError:  # pragma: no cover - build-dependent
    _rdf_kernel = None

rdf_module = sys.modules[RDF.__module__]

WATER = np.array([[0.0, 0, 0], [0.9572, 0, 0], [-0.2399, 0.9266, 0]])
METHANE = np.array([
    [0.0, 0, 0],
    [0.63, 0.63, 0.63],
    [-0.63, -0.63, 0.63],
    [-0.63, 0.63, -0.63],
    [0.63, -0.63, -0.63],
])

#: (name, element symbols, geometry, number of molecules)
MOLECULES = [
    ("WAT", ("O", "H", "H"), WATER, 24),
    ("MET", ("C", "H", "H", "H", "H"), METHANE, 6),
    ("ION", ("Na",), np.zeros((1, 3)), 4),
]

N_FRAMES = 4
DELTA_R = 0.1
R_MAX = 4.5
N_BINS = 45
BOX = 10.0

#: selections of the atoms of the three kinds of molecules
WATERS = "0-71"
METHANES = "72-101"
IONS = "102-105"


def _masses(symbols):
    return np.array([Atom(symbol).mass for symbol in symbols])


def build_topology():
    """Waters (atoms 0-71), methanes (72-101) and ions (102-105)."""

    residues, atoms, residue_ids = [], [], []
    residue_id = 0

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
            residue_id += 1
            atoms += [Atom(symbol) for symbol in symbols]
            residue_ids += [index + 1] * len(symbols)

    return Topology(
        atoms=atoms,
        residue_ids=np.array(residue_ids),
        reference_residues=residues,
    )


def molecule_ids(topology):
    return topology.residue_numbers


def generate_frames(box_sizes, seed=7):
    """
    Seeded frames with wrapped atoms; values are rounded to the digits of the
    xyz files, so that files and arrays hold exactly the same numbers.
    """

    rng = np.random.default_rng(seed)
    frames = []

    for box in box_sizes:
        molecules = []

        for _, symbols, geometry, count in MOLECULES:
            masses = _masses(symbols)
            geometry = geometry - (masses[:, None] * geometry).sum(0) / masses.sum()

            for _ in range(count):
                rotation = Rotation.from_quat(rng.normal(size=4))
                center = rng.uniform(0.0, box, 3)
                molecules.append(rotation.apply(geometry) + center)

        positions = np.vstack(molecules) % box
        frames.append(
            np.array([[float(f"{value:.6f}") for value in row] for row in positions])
        )

    return frames


def write_xyz(path, frames, box_sizes, topology):
    names = [atom.name for atom in topology.atoms]
    lines = []

    for positions, box in zip(frames, box_sizes):
        lines.append(f"{len(names)} {box} {box} {box} 90.0 90.0 90.0\n\n")
        lines += [
            f"{name} {x:.6f} {y:.6f} {z:.6f}\n"
            for name, (x, y, z) in zip(names, positions)
        ]

    path.write_text("".join(lines), encoding="utf-8")

    return str(path)


def selected_atoms(spec):
    """Atom indices of a selection written as comma separated ranges."""

    atoms = []
    for part in spec.split(","):
        first, _, last = part.partition("-")
        atoms += list(range(int(first), int(last or first) + 1))

    return np.array(atoms)


def centers_of_mass(positions, atoms, molecule, masses, box):
    """True centers of mass of the molecules of the given atoms."""

    ids = np.unique(molecule[atoms])
    coms = []

    for molecule_id in ids:
        chosen = atoms[molecule[atoms] == molecule_id]
        relative = positions[chosen] - positions[chosen[0]]
        relative -= box * np.round(relative / box)
        coms.append(
            positions[chosen[0]]
            + (relative * masses[chosen, None]).sum(0) / masses[chosen].sum()
        )

    return ids, np.array(coms)


def brute_force_histogram(frames, box_sizes, topology, reference, target):
    """
    Histogram of the minimum image distances between the centers of mass
    of the reference and the target molecules, without the pair of a
    molecule with itself.
    """

    molecule = molecule_ids(topology)
    masses = np.array([atom.mass for atom in topology.atoms])
    histogram = np.zeros(N_BINS)

    for positions, box in zip(frames, box_sizes):
        ref_ids, ref = centers_of_mass(
            positions, selected_atoms(reference), molecule, masses, box
        )
        tgt_ids, tgt = centers_of_mass(
            positions, selected_atoms(target), molecule, masses, box
        )

        difference = ref[:, None] - tgt[None]
        difference -= box * np.round(difference / box)
        distance = np.linalg.norm(difference, axis=-1)
        distance[ref_ids[:, None] == tgt_ids[None]] = np.inf

        bins = np.floor(distance / DELTA_R)
        bins = bins[bins < N_BINS].astype(int)
        histogram += np.bincount(bins, minlength=N_BINS)

    return histogram


def in_memory(frames, box_sizes, topology):
    return Trajectory([
        AtomicSystem(
            topology=topology,
            pos=positions,
            cell=Cell(box, box, box),
        ) for positions, box in zip(frames, box_sizes)
    ])


def run_com_rdf(traj, reference, target, **kwargs):
    analysis = RDF(
        traj,
        reference,
        target,
        delta_r=DELTA_R,
        r_max=R_MAX,
        reference_com=True,
        target_com=True,
        **kwargs,
    )
    analysis.run()

    return analysis


BOXES = [BOX] * N_FRAMES
NPT_BOXES = [BOX + 0.1 * index for index in range(N_FRAMES)]

#: (reference, target): all molecules, disjoint kinds, shared molecules,
#: and the hydrogens only (centers over a part of the atoms)
SELECTIONS = {
    "all": ("0-105", "0-105"),
    "ions-around-waters": (IONS, WATERS),
    "waters-around-waters-and-methanes": (WATERS, "0-101"),
    "waters-and-methanes-with-the-same-atoms": ("0-101", "0-101"),
}


@pytest.fixture(name="topology")
def topology_fixture():
    return build_topology()


@pytest.mark.parametrize("boxes", [BOXES, NPT_BOXES], ids=["nvt", "npt"])
@pytest.mark.parametrize(
    "selection", SELECTIONS.values(), ids=SELECTIONS.keys()
)
class TestAgainstBruteForce:

    def test_in_memory_trajectory(self, topology, selection, boxes):
        frames = generate_frames(boxes)
        reference, target = selection

        analysis = run_com_rdf(
            in_memory(frames, boxes, topology), reference, target
        )

        assert analysis._raw_reader is None
        assert np.array_equal(
            analysis.bins,
            brute_force_histogram(frames, boxes, topology, reference, target),
        )
        assert analysis.bins.sum() > 0


FILE_PATHS = [
    pytest.param(
        _rdf_kernel_py.rdf_frame_histogram, "stream", id="fallback-stream"
    ),
    pytest.param(
        _rdf_kernel.rdf_frame_histogram if _rdf_kernel else None,
        "stream",
        id="cython-stream",
        marks=pytest.mark.skipif(
            _rdf_kernel is None, reason="Cython _rdf_kernel not built"
        ),
    ),
    pytest.param(
        _rdf_kernel.rdf_frame_histogram if _rdf_kernel else None,
        "batch",
        id="cython-batch",
        marks=pytest.mark.skipif(
            _rdf_kernel is None, reason="Cython _rdf_kernel not built"
        ),
    ),
    pytest.param(
        _rdf_kernel.rdf_frame_histogram if _rdf_kernel else None,
        "parallel",
        id="cython-batch-parallel",
        marks=pytest.mark.skipif(
            _rdf_kernel is None, reason="Cython _rdf_kernel not built"
        ),
    ),
]


@pytest.mark.parametrize("boxes", [BOXES, NPT_BOXES], ids=["nvt", "npt"])
@pytest.mark.parametrize("kernel, path", FILE_PATHS)
class TestFilePaths:

    @pytest.fixture(autouse=True)
    def configure(self, kernel, path, monkeypatch):
        monkeypatch.setattr(rdf_module, "rdf_frame_histogram", kernel)

        if path == "stream":
            monkeypatch.setattr(RDF, "_batch_max_bytes", 0)
        if path == "parallel":
            monkeypatch.setattr(RDF, "_batch_pairs_per_worker", 1)
            monkeypatch.setattr(RDF, "_batch_max_workers", 2)
            monkeypatch.setattr(rdf_module, "cpu_count", lambda: 2)

        self.batch_calls = []
        histogram_range = rdf_module._raw_histogram_range

        def spy(*arguments):
            self.batch_calls.append(1)
            return histogram_range(*arguments)

        monkeypatch.setattr(rdf_module, "_raw_histogram_range", spy)

    @pytest.mark.parametrize(
        "selection", SELECTIONS.values(), ids=SELECTIONS.keys()
    )
    def test_matches_brute_force(
        self, topology, selection, boxes, path, tmp_path
    ):
        frames = generate_frames(boxes)
        reference, target = selection
        filename = write_xyz(tmp_path / "traj.xyz", frames, boxes, topology)

        analysis = run_com_rdf(
            TrajectoryReader(filename, topology=topology), reference, target
        )

        assert analysis._raw_reader is not None
        assert len(self.batch_calls) == {"stream": 0, "batch": 1, "parallel": 2}[path]
        assert np.array_equal(
            analysis.bins,
            brute_force_histogram(frames, boxes, topology, reference, target),
        )

    def test_equals_the_in_memory_trajectory(
        self, topology, boxes, tmp_path
    ):
        frames = generate_frames(boxes)
        filename = write_xyz(tmp_path / "traj.xyz", frames, boxes, topology)

        from_file = run_com_rdf(
            TrajectoryReader(filename, topology=topology), "0-105", "0-105"
        )
        from_memory = run_com_rdf(
            in_memory(frames, boxes, topology), "0-105", "0-105"
        )

        assert np.array_equal(from_file.bins, from_memory.bins)
        assert np.array_equal(
            from_file.normalized_bins, from_memory.normalized_bins
        )


class TestNormalization:

    """
    The reference count and the target density are numbers of molecules, the
    coordination number counts molecules.
    """

    N_MOLECULES = 200
    BOX = 12.0
    FRAMES = 30

    @pytest.fixture(name="ideal_gas")
    def ideal_gas_fixture(self):
        residue = Residue(
            name="WAT",
            elements=[Element("O"), Element("H"), Element("H")],
            residue_id=1,
            total_charge=0,
            atom_types=np.array([0, 1, 1]),
            partial_charges=np.zeros(3),
        )
        topology = Topology(
            atoms=[Atom(symbol) for symbol in "OHH" * self.N_MOLECULES],
            residue_ids=np.ones(3 * self.N_MOLECULES, dtype=int),
            reference_residues=[residue],
        )
        masses = _masses("OHH")
        geometry = WATER - (masses[:, None] * WATER).sum(0) / masses.sum()
        rng = np.random.default_rng(11)
        cell = Cell(self.BOX, self.BOX, self.BOX)

        frames = []
        for _ in range(self.FRAMES):
            centers = rng.uniform(0.0, self.BOX, (self.N_MOLECULES, 1, 3))
            rotations = Rotation.from_quat(rng.normal(size=(self.N_MOLECULES, 4)))
            molecules = np.stack([rotations[i].apply(geometry) for i in range(self.N_MOLECULES)])
            frames.append(
                AtomicSystem(
                    topology=topology,
                    pos=((molecules + centers).reshape(-1, 3)) % self.BOX,
                    cell=cell,
                )
            )

        return Trajectory(frames)

    def test_the_rdf_of_uniformly_placed_molecules_is_one(self, ideal_gas):
        analysis = RDF(
            ideal_gas, "*", "*", delta_r=0.5, r_max=5.5,
            reference_com=True, target_com=True,
        )
        distances, g, coordination, *_ = analysis.run()

        n = self.N_MOLECULES
        # the pair of a molecule with itself is excluded, the density is
        # N / V as for atoms, so the plateau is (N - 1) / N
        plateau = g[distances > 1.5]
        assert np.allclose(plateau, (n - 1) / n, atol=0.04)
        assert abs(plateau.mean() - (n - 1) / n) < 0.01

        volume = self.BOX**3
        expected = (n - 1) * 4.0 / 3.0 * np.pi * 5.5**3 / volume
        assert coordination[-1] == pytest.approx(expected, rel=0.01)

        # normalizing by atoms would be wrong by a factor of three
        assert len(analysis.reference_indices) == n
        assert len(analysis.target_indices) == n


class TestOptions:

    @staticmethod
    def traj(topology):
        return in_memory(generate_frames([BOX] * 2), [BOX] * 2, topology)

    def test_no_intra_molecular_has_no_effect(self, topology):
        plain = run_com_rdf(self.traj(topology), "0-105", "0-105")
        excluded = run_com_rdf(
            self.traj(topology), "0-105", "0-105", no_intra_molecular=True
        )

        assert excluded.no_intra_molecular is False
        assert np.array_equal(plain.bins, excluded.bins)
        assert np.array_equal(plain.normalized_bins, excluded.normalized_bins)

    def test_the_file_fast_path_is_kept_with_no_intra_molecular(
        self, topology, tmp_path
    ):
        frames = generate_frames([BOX] * 2)
        filename = write_xyz(tmp_path / "traj.xyz", frames, [BOX] * 2, topology)

        analysis = run_com_rdf(
            TrajectoryReader(filename, topology=topology),
            "0-105",
            "0-105",
            no_intra_molecular=True,
        )

        assert analysis._raw_reader is not None

    def test_the_reference_indices_are_the_rows_of_the_centers(self, topology):
        analysis = RDF(
            self.traj(topology), IONS, WATERS, delta_r=DELTA_R, r_max=R_MAX,
            reference_com=True, target_com=True,
        )

        n_atoms = topology.n_atoms
        assert len(analysis.reference_indices) == 4
        assert len(analysis.target_indices) == 24
        assert analysis.reference_indices.min() >= n_atoms
        assert analysis.target_indices.min() >= n_atoms
        assert not set(analysis.reference_indices) & set(analysis.target_indices)

    def test_shared_molecules_share_their_rows(self, topology):
        analysis = RDF(
            self.traj(topology), WATERS, "0-101", delta_r=DELTA_R, r_max=R_MAX,
            reference_com=True, target_com=True,
        )

        assert set(analysis.reference_indices) <= set(analysis.target_indices)
        assert len(analysis.target_indices) == 30

    def test_a_molecule_without_selected_atoms_is_left_out(self, topology):
        analysis = RDF(
            self.traj(topology), IONS, "Na", delta_r=DELTA_R, r_max=R_MAX,
            reference_com=True, target_com=True,
        )

        assert len(analysis.reference_indices) == 4
        assert len(analysis.target_indices) == 4


class TestErrors:

    @staticmethod
    def traj(topology):
        return in_memory(generate_frames([BOX] * 2), [BOX] * 2, topology)

    def test_a_topology_without_residues_raises(self):
        atoms = [Atom(symbol) for symbol in "OHH" * 4]
        frames = [
            AtomicSystem(
                atoms=atoms,
                pos=generate_frames([BOX])[0][:12],
                cell=Cell(BOX, BOX, BOX),
            )
        ]

        with pytest.raises(RDFError, match="topology has no residues"):
            RDF(
                Trajectory(frames), "O", "H", delta_r=DELTA_R, r_max=R_MAX,
                reference_com=True, target_com=True,
            )

    def test_a_shared_molecule_with_different_atoms_raises(self, topology):
        # the oxygens define the reference center, all atoms the target one
        with pytest.raises(RDFError, match="needs the same selected atoms"):
            RDF(
                self.traj(topology), "O", WATERS, delta_r=DELTA_R, r_max=R_MAX,
                reference_com=True, target_com=True,
            )

    def test_atoms_without_a_mass_raise(self):
        residue = Residue(
            name="MOL",
            elements=[Element("C"), Element("C")],
            residue_id=1,
            total_charge=0,
            atom_types=np.array([0, 0]),
            partial_charges=np.zeros(2),
        )
        topology = Topology(
            atoms=[Atom(name, use_guess_element=False) for name in ("A1", "A2")] * 2,
            residue_ids=np.ones(4, dtype=int),
            reference_residues=[residue],
        )
        traj = Trajectory([
            AtomicSystem(
                topology=topology,
                pos=np.arange(12, dtype=float).reshape(4, 3),
                cell=Cell(BOX, BOX, BOX),
            )
        ])

        with pytest.raises(RDFError, match="has no mass|no mass"):
            RDF(
                traj, "*", "*", delta_r=DELTA_R, r_max=R_MAX,
                reference_com=True, target_com=True,
            )

    def test_a_molecule_larger_than_half_the_box_raises(self):
        residue = Residue(
            name="LONG",
            elements=[Element("C")] * 3,
            residue_id=1,
            total_charge=0,
            atom_types=np.zeros(3, dtype=int),
            partial_charges=np.zeros(3),
        )
        topology = Topology(
            atoms=[Atom("C")] * 6,
            residue_ids=np.ones(6, dtype=int),
            reference_residues=[residue],
        )
        positions = np.array([
            [0.0, 0, 0], [3.2, 0, 0], [6.4, 0, 0],
            [0.0, 5, 0], [0.5, 5, 0], [1.0, 5, 0],
        ])
        traj = Trajectory([
            AtomicSystem(
                topology=topology, pos=positions, cell=Cell(BOX, BOX, BOX)
            )
        ])
        analysis = RDF(
            traj, "*", "*", delta_r=DELTA_R, r_max=R_MAX,
            reference_com=True, target_com=True,
        )

        with pytest.raises(AtomicSystemError, match="spans half of the box"):
            analysis.run()

    def test_a_file_frame_with_other_atoms_than_the_topology_raises(
        self, topology, tmp_path, monkeypatch
    ):
        frames = generate_frames([BOX] * 2)
        filename = write_xyz(tmp_path / "traj.xyz", frames, [BOX] * 2, topology)
        analysis = RDF(
            TrajectoryReader(filename, topology=topology), WATERS, WATERS,
            delta_r=DELTA_R, r_max=R_MAX,
            reference_com=True, target_com=True,
        )

        # the waters (atoms 0-71) are all there, but the frame has 100 atoms
        def truncated_generator():
            yield frames[0][:100], Cell(BOX, BOX, BOX)

        monkeypatch.setattr(
            analysis._raw_reader, "raw_frame_generator", truncated_generator
        )
        monkeypatch.setattr(analysis, "_batch_max_bytes", 0)

        with pytest.raises(RDFError, match="has 100 atoms, but the topology"):
            analysis.run()


class TestPlaceholderAtoms:

    def test_a_placeholder_mass_logs_a_warning(self, caplog):
        residue = Residue(
            name="MOL",
            elements=[Element("O"), Element("H"), Element("X")],
            residue_id=1,
            total_charge=0,
            atom_types=np.arange(3),
            partial_charges=np.zeros(3),
        )
        topology = Topology(
            atoms=[Atom(name) for name in "OHX" * 4],
            residue_ids=np.ones(12, dtype=int),
            reference_residues=[residue],
        )
        traj = Trajectory([
            AtomicSystem(
                topology=topology,
                pos=generate_frames([BOX])[0][:12],
                cell=Cell(BOX, BOX, BOX),
            )
        ])

        with caplog.at_level("WARNING"):
            RDF(
                traj, "*", "*", delta_r=DELTA_R, r_max=R_MAX,
                reference_com=True, target_com=True,
            )

        assert "placeholder mass (x)" in caplog.text.lower()


class TestLog:

    def test_the_log_counts_molecules(self, topology, tmp_path):
        from PQAnalysis.analysis.rdf import RDFLogWriter

        analysis = RDF(
            in_memory(generate_frames([BOX]), [BOX], topology),
            IONS, WATERS, delta_r=DELTA_R, r_max=R_MAX,
            reference_com=True, target_com=True,
        )

        log = tmp_path / "rdf.log"
        writer = RDFLogWriter(str(log))
        writer.write_before_run(analysis)

        text = log.read_text(encoding="utf-8")
        assert "total number of molecules (centers of mass) in reference selection: 4" in text
        assert "total number of molecules (centers of mass) in target selection:    24" in text



def brute_force_mixed(
    frames, boxes, topology, reference, target, flags, exclude_intra
):
    """
    Histogram and mean number of targets per reference for the RDF between
    centers and atoms. A side with a center uses the true centers of mass of
    the molecules of its selected atoms, the other side its selected atoms.
    """

    reference_com, target_com = flags
    molecule = molecule_ids(topology)
    masses = np.array([atom.mass for atom in topology.atoms])
    histogram = np.zeros(N_BINS)
    n_targets = []

    def entities(positions, spec, as_center, box):
        atoms = selected_atoms(spec)
        if as_center:
            return centers_of_mass(positions, atoms, molecule, masses, box)
        return molecule[atoms], positions[atoms]

    for positions, box in zip(frames, boxes):
        ref_molecules, ref = entities(positions, reference, reference_com, box)
        tgt_molecules, tgt = entities(positions, target, target_com, box)

        difference = ref[:, None] - tgt[None]
        difference -= box * np.round(difference / box)
        distance = np.linalg.norm(difference, axis=-1)
        same = ref_molecules[:, None] == tgt_molecules[None]
        if exclude_intra:
            distance[same] = np.inf
            n_targets.append((~same).sum(axis=1).mean())
        else:
            n_targets.append(len(tgt))

        bins = np.floor(distance / DELTA_R)
        bins = bins[bins < N_BINS].astype(int)
        histogram += np.bincount(bins, minlength=N_BINS)

    return histogram, np.mean(n_targets)


#: (reference, target, reference_com, target_com): atoms around centers,
#: centers around atoms, atoms of two kinds around centers of a third, and
#: centers of molecules of different size around all atoms
MIXED = {
    "oxygens-around-water-centers": ("O", "0-71", False, True),
    "ions-around-all-centers": (IONS, "0-105", False, True),
    "water-centers-around-hydrogens": (WATERS, "H", True, False),
    "all-centers-around-all-atoms": ("0-105", "0-105", True, False),
    # centers of a subset that does not start with the first molecule: the
    # row of a center is not its molecule number
    "methane-centers-around-all-atoms": (METHANES, "0-105", True, False),
    "all-atoms-around-ion-centers": ("0-105", IONS, False, True),
}


def run_mixed(traj, reference, target, flags, **kwargs):
    analysis = RDF(
        traj,
        reference,
        target,
        delta_r=DELTA_R,
        r_max=R_MAX,
        reference_com=flags[0],
        target_com=flags[1],
        **kwargs,
    )
    analysis.run()

    return analysis


def mixed_spec(spec):
    """Selection string for the RDF and for the brute force."""

    if spec == "O":
        return "O", ",".join(
            f"{3 * i}" for i in range(24)
        )
    if spec == "H":
        # the hydrogens of the waters and of the methanes
        return "H", ",".join(
            [f"{3 * i + 1}-{3 * i + 2}" for i in range(24)]
            + [f"{72 + 5 * i + 1}-{72 + 5 * i + 4}" for i in range(6)]
        )
    return spec, spec


@pytest.mark.parametrize("exclude_intra", [False, True], ids=["intra", "no-intra"])
@pytest.mark.parametrize("boxes", [BOXES, NPT_BOXES], ids=["nvt", "npt"])
@pytest.mark.parametrize("case", MIXED.values(), ids=MIXED.keys())
class TestAtomsAndCenters:

    def test_in_memory_matches_brute_force(
        self, topology, case, boxes, exclude_intra
    ):
        reference, target, *flags = case
        frames = generate_frames(boxes)
        rdf_ref, brute_ref = mixed_spec(reference)
        rdf_tgt, brute_tgt = mixed_spec(target)

        analysis = run_mixed(
            in_memory(frames, boxes, topology), rdf_ref, rdf_tgt, flags,
            no_intra_molecular=exclude_intra,
        )
        histogram, mean_targets = brute_force_mixed(
            frames, boxes, topology, brute_ref, brute_tgt, flags,
            exclude_intra,
        )

        assert np.array_equal(analysis.bins, histogram)
        assert histogram.sum() > 0
        if exclude_intra:
            assert analysis._raw_reader is None
            expected = (
                np.mean([len(t) for t in analysis.target_index_combinations])
            )
            assert expected == pytest.approx(mean_targets)

    def test_the_file_matches_the_in_memory_trajectory(
        self, topology, case, boxes, exclude_intra, tmp_path
    ):
        reference, target, *flags = case
        frames = generate_frames(boxes)
        rdf_ref, _ = mixed_spec(reference)
        rdf_tgt, _ = mixed_spec(target)
        filename = write_xyz(tmp_path / "traj.xyz", frames, boxes, topology)

        from_file = run_mixed(
            TrajectoryReader(filename, topology=topology), rdf_ref, rdf_tgt,
            flags, no_intra_molecular=exclude_intra,
        )
        from_memory = run_mixed(
            in_memory(frames, boxes, topology), rdf_ref, rdf_tgt, flags,
            no_intra_molecular=exclude_intra,
        )

        # without the exclusion the file takes the kernel fast path
        assert (from_file._raw_reader is not None) == (not exclude_intra)
        assert np.array_equal(from_file.bins, from_memory.bins)
        assert np.array_equal(
            from_file.normalized_bins, from_memory.normalized_bins
        )


@pytest.mark.parametrize("kernel, path", FILE_PATHS)
class TestAtomsAndCentersFilePaths:

    @pytest.fixture(autouse=True)
    def configure(self, kernel, path, monkeypatch):
        monkeypatch.setattr(rdf_module, "rdf_frame_histogram", kernel)

        if path == "stream":
            monkeypatch.setattr(RDF, "_batch_max_bytes", 0)
        if path == "parallel":
            monkeypatch.setattr(RDF, "_batch_pairs_per_worker", 1)
            monkeypatch.setattr(RDF, "_batch_max_workers", 2)
            monkeypatch.setattr(rdf_module, "cpu_count", lambda: 2)

    @pytest.mark.parametrize("case", MIXED.values(), ids=MIXED.keys())
    def test_matches_brute_force(self, topology, case, path, tmp_path):
        reference, target, *flags = case
        frames = generate_frames(NPT_BOXES)
        rdf_ref, brute_ref = mixed_spec(reference)
        rdf_tgt, brute_tgt = mixed_spec(target)
        filename = write_xyz(
            tmp_path / "traj.xyz", frames, NPT_BOXES, topology
        )

        analysis = run_mixed(
            TrajectoryReader(filename, topology=topology), rdf_ref, rdf_tgt,
            flags,
        )
        histogram, _ = brute_force_mixed(
            frames, NPT_BOXES, topology, brute_ref, brute_tgt, flags, False
        )

        assert analysis._raw_reader is not None
        assert np.array_equal(analysis.bins, histogram)


class TestAtomsAndCentersNormalization:

    """
    Oxygens around the centers of 200 uniformly placed waters: the oxygen of
    a molecule is 0.065 A from its own center, every other center is
    uncorrelated with it.
    """

    N_MOLECULES = TestNormalization.N_MOLECULES
    BOX = TestNormalization.BOX
    FRAMES = TestNormalization.FRAMES

    @pytest.fixture(name="ideal_gas")
    def ideal_gas_fixture(self):
        return TestNormalization().ideal_gas_fixture.__wrapped__(self)

    def run(self, ideal_gas, no_intra_molecular):
        return RDF(
            ideal_gas, "O", "*", delta_r=0.5, r_max=5.5,
            reference_com=False, target_com=True,
            no_intra_molecular=no_intra_molecular,
        )

    def test_the_own_center_is_excluded_from_the_density(self, ideal_gas):
        analysis = self.run(ideal_gas, True)
        distances, g, coordination, *_ = analysis.run()

        n = self.N_MOLECULES
        # every oxygen has N - 1 centers to look at: the density is the
        # density of the others, so the plateau is exactly one
        assert all(len(t) == n - 1 for t in analysis.target_index_combinations)
        # the first bin only holds the few other centers that are there by
        # chance (about 360 of 6000 pairs), not the own centers
        assert analysis.bins[0] < 0.1 * n * self.FRAMES
        assert g[0] == pytest.approx(1.0, abs=0.25)
        assert np.allclose(g[distances > 1.5], 1.0, atol=0.04)
        expected = (n - 1) * 4.0 / 3.0 * np.pi * 5.5**3 / self.BOX**3
        assert coordination[-1] == pytest.approx(expected, rel=0.01)

    def test_the_own_center_is_a_pair_without_the_exclusion(self, ideal_gas):
        analysis = self.run(ideal_gas, False)
        distances, g, coordination, *_ = analysis.run()

        n = self.N_MOLECULES
        # the own center is 0.065 A away: it fills the first bin and the
        # plateau is that of N - 1 targets in a density of N
        own = self.N_MOLECULES * self.FRAMES
        assert own < analysis.bins[0] < 1.1 * own
        assert np.allclose(g[distances > 1.5], (n - 1) / n, atol=0.04)
        assert coordination[0] == pytest.approx(1.0, abs=0.1)

    def test_the_mean_of_the_excluded_atoms_is_the_density(self, topology):
        # centers of waters (3 atoms) and methanes (5 atoms) around all
        # hydrogens: a water center excludes 2 of them, a methane center 4
        boxes = [BOX] * 2
        analysis = run_mixed(
            in_memory(generate_frames(boxes), boxes, topology),
            "0-101", "H", (True, False), no_intra_molecular=True,
        )

        counts = {len(t) for t in analysis.target_index_combinations}
        n_hydrogens = 24 * 2 + 6 * 4
        assert counts == {n_hydrogens - 2, n_hydrogens - 4}

        _, brute_mean = brute_force_mixed(
            generate_frames(boxes), boxes, topology, "0-101",
            ",".join(f"{3 * i + 1}-{3 * i + 2}" for i in range(24))
            + "," + ",".join(f"{72 + 5 * i + 1}-{72 + 5 * i + 4}" for i in range(6)),
            (True, False), True,
        )
        mean = np.mean([len(t) for t in analysis.target_index_combinations])
        assert mean == pytest.approx(brute_mean)
        density = mean / analysis.average_volume
        norm = analysis._norm(
            analysis.n_bins, analysis.delta_r, density,
            len(analysis.reference_indices), analysis.n_frames, 0.0,
        )
        assert np.allclose(analysis.normalized_bins, analysis.bins / norm)


class TestAtomsAndCentersOptions:

    @staticmethod
    def traj(topology):
        return in_memory(generate_frames([BOX] * 2), [BOX] * 2, topology)

    def test_the_side_without_centers_keeps_its_atoms(self, topology):
        analysis = RDF(
            self.traj(topology), IONS, WATERS, delta_r=DELTA_R, r_max=R_MAX,
            reference_com=False, target_com=True,
        )

        assert list(analysis.reference_indices) == [102, 103, 104, 105]
        assert analysis.target_indices.min() >= topology.n_atoms
        assert len(analysis.target_indices) == 24

    def test_only_the_atoms_of_the_center_side_need_a_mass(self):
        residue = Residue(
            name="MOL",
            elements=[Element("C"), Element("C")],
            residue_id=1,
            total_charge=0,
            atom_types=np.array([0, 0]),
            partial_charges=np.zeros(2),
        )
        atoms = [Atom("C"), Atom("A1", use_guess_element=False)] * 2
        topology = Topology(
            atoms=atoms, residue_ids=np.ones(4, dtype=int),
            reference_residues=[residue],
        )
        traj = Trajectory([
            AtomicSystem(
                topology=topology,
                pos=np.arange(12, dtype=float).reshape(4, 3),
                cell=Cell(BOX, BOX, BOX),
            )
        ])

        # the massless atoms are only targets, the centers use the carbons
        RDF(
            traj, "C", "*", delta_r=DELTA_R, r_max=R_MAX,
            reference_com=True, target_com=False,
        )
        with pytest.raises(RDFError, match="no mass"):
            RDF(
                traj, "*", "C", delta_r=DELTA_R, r_max=R_MAX,
                reference_com=True, target_com=False,
            )

    def test_a_file_frame_that_lacks_the_atoms_of_the_atom_side_raises(
        self, topology, tmp_path, monkeypatch
    ):
        frames = generate_frames([BOX] * 2)
        filename = write_xyz(tmp_path / "traj.xyz", frames, [BOX] * 2, topology)
        analysis = RDF(
            TrajectoryReader(filename, topology=topology), IONS, WATERS,
            delta_r=DELTA_R, r_max=R_MAX,
            reference_com=False, target_com=True,
        )

        # the waters that define the centers are there (atoms 0-71), the
        # ions that are the references (atoms 102-105) are not
        def truncated_generator():
            yield frames[0][:100], Cell(BOX, BOX, BOX)

        monkeypatch.setattr(
            analysis._raw_reader, "raw_frame_generator", truncated_generator
        )
        monkeypatch.setattr(analysis, "_batch_max_bytes", 0)

        with pytest.raises(RDFError, match="provides only 100 atoms"):
            analysis.run()

    def test_atoms_of_the_same_molecule_can_differ_from_the_center_atoms(
        self, topology
    ):
        # no consistency rule for one-sided use: a water oxygen and the
        # center of the same water are different entities
        analysis = run_mixed(
            self.traj(topology), "O", WATERS, (False, True),
            no_intra_molecular=False,
        )

        assert analysis.bins.sum() > 0

    def test_a_molecule_larger_than_half_the_box_raises_when_running(self):
        residue = Residue(
            name="LONG",
            elements=[Element("C")] * 3,
            residue_id=1,
            total_charge=0,
            atom_types=np.zeros(3, dtype=int),
            partial_charges=np.zeros(3),
        )
        topology = Topology(
            atoms=[Atom("C")] * 3, residue_ids=np.ones(3, dtype=int),
            reference_residues=[residue],
        )
        traj = Trajectory([
            AtomicSystem(
                topology=topology,
                pos=np.array([[0.0, 0, 0], [3.2, 0, 0], [6.4, 0, 0]]),
                cell=Cell(BOX, BOX, BOX),
            )
        ])
        analysis = RDF(
            traj, "*", "*", delta_r=DELTA_R, r_max=R_MAX,
            reference_com=False, target_com=True,
        )

        with pytest.raises(AtomicSystemError, match="spans half of the box"):
            analysis.run()

    def test_the_log_names_the_entities_of_each_side(self, topology, tmp_path):
        from PQAnalysis.analysis.rdf import RDFLogWriter

        analysis = RDF(
            self.traj(topology), IONS, WATERS, delta_r=DELTA_R, r_max=R_MAX,
            reference_com=False, target_com=True,
        )
        log = tmp_path / "rdf.log"
        RDFLogWriter(str(log)).write_before_run(analysis)

        text = log.read_text(encoding="utf-8")
        assert "total number of atoms in reference selection: 4" in text
        assert "total number of molecules (centers of mass) in target selection:    24" in text
