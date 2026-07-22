"""
Equivalence tests of the Cython RDF distance-histogram kernel against
the pure Python/numpy fallback kernel and bit-identity tests of the
raw-frame fast path against the original RDF path.

Both kernels implement the identical signature and are driven frame
by frame over synthetic float32 frames (orthorhombic, triclinic, NPT
with changing boxes, vacuum and bin-edge cases). The resulting int64
histograms must agree exactly (atol 0). The fast path of the RDF
class itself is additionally run against both kernels via
monkeypatching and compared bit-identically (np.array_equal on all
output arrays) against the original path, both from an in-memory
trajectory and from a TrajectoryReader with the fast path disabled.
"""

import sys

import numpy as np
import pytest

from PQAnalysis.analysis import RDF
from PQAnalysis.analysis.rdf import _rdf_kernel_py
from PQAnalysis.analysis.rdf.exceptions import RDFError
from PQAnalysis.core import Cell
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.io.traj_file.exceptions import TrajectoryReaderError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception

# pylint: disable=protected-access

try:
    from PQAnalysis.analysis.rdf import _rdf_kernel
except ModuleNotFoundError:  # pragma: no cover - build-dependent
    _rdf_kernel = None

#: The module defining the RDF class (the package attribute ``rdf``
#: is shadowed by the api function of the same name).
rdf_module = sys.modules[RDF.__module__]

KERNELS = [
    pytest.param(_rdf_kernel_py.rdf_frame_histogram, id="python-fallback"),
    pytest.param(
        _rdf_kernel.rdf_frame_histogram if _rdf_kernel is not None else None,
        id="cython",
        marks=pytest.mark.skipif(
            _rdf_kernel is None,
            reason="Cython _rdf_kernel extension not built",
        ),
    ),
]


def _random_frames(n_frames, n_atoms, spread, seed):
    """
    Builds float32 frames of uniformly distributed positions.
    """
    rng = np.random.default_rng(seed)

    return [
        np.asarray(
            rng.uniform(0.0, spread, size=(n_atoms, 3)),
            dtype=np.float32,
        )
        for _ in range(n_frames)
    ]


def _drive_kernel(kernel, frames, cells, r_min, delta_r, n_bins):
    """
    Drives a kernel implementation frame by frame (like the fast
    path of the RDF class does) and returns the histogram.
    """
    reference_indices = np.arange(0, 30, 2, dtype=np.int64)
    # overlaps with the reference indices, so that the self-pair
    # exclusion is exercised
    target_indices = np.arange(0, 30, 3, dtype=np.int64)

    hist = np.zeros(n_bins, dtype=np.int64)

    for values, cell in zip(frames, cells):
        is_orthorhombic = 1 if (
            cell.alpha == 90 and cell.beta == 90 and cell.gamma == 90
        ) else 0

        kernel(
            values,
            reference_indices,
            target_indices,
            np.ascontiguousarray(cell.box_lengths, dtype=np.float64),
            np.ascontiguousarray(cell.box_matrix, dtype=np.float64),
            np.ascontiguousarray(cell.inverse_box_matrix, dtype=np.float64),
            is_orthorhombic,
            r_min,
            delta_r,
            n_bins,
            hist,
        )

    return hist


def _assert_kernel_equivalence(cells, r_min=0.0, delta_r=0.05, n_bins=200):
    """
    Runs both kernels over the same frames and asserts that the
    histograms agree exactly.
    """
    if _rdf_kernel is None:  # pragma: no cover - build-dependent
        pytest.skip("Cython _rdf_kernel extension not built")

    frames = _random_frames(len(cells), 30, 17.0, seed=2026)

    histograms = [
        _drive_kernel(kernel, frames, cells, r_min, delta_r, n_bins)
        for kernel in (
            _rdf_kernel.rdf_frame_histogram,
            _rdf_kernel_py.rdf_frame_histogram,
        )
    ]

    assert np.array_equal(histograms[0], histograms[1])

    # both kernels must have accumulated a non-trivial histogram
    assert histograms[0].sum() > 0


def _write_trajectory(path, positions, headers):
    """
    Writes an xyz trajectory file with the given per-frame header
    box strings.
    """
    with open(path, "w", encoding="utf-8") as file:
        for frame_positions, header in zip(positions, headers):
            file.write(f"{len(frame_positions)}{header}\n\n")

            for i, (x, y, z) in enumerate(frame_positions):
                name = "O" if i % 2 == 0 else "H"
                file.write(f"{name} {x:.6f} {y:.6f} {z:.6f}\n")

    return str(path)


def _write_random_trajectory(path, n_frames=25, n_atoms=16, seed=4242):
    """
    Writes a constant orthorhombic box xyz trajectory file.
    """
    rng = np.random.default_rng(seed)

    positions = rng.uniform(0.0, 12.0, size=(n_frames, n_atoms, 3))
    headers = [" 12.0 12.0 12.0"] * n_frames

    return _write_trajectory(path, positions, headers)


def _write_npt_trajectory(path, n_frames=24, n_atoms=16, seed=99):
    """
    Writes an NPT-like xyz trajectory file with per-frame changing
    boxes, mixing orthorhombic and triclinic headers.
    """
    rng = np.random.default_rng(seed)

    positions = rng.uniform(0.0, 12.0, size=(n_frames, n_atoms, 3))

    headers = []
    for i in range(n_frames):
        factor = 1.0 + 0.01 * (i % 4)
        if (i // 4) % 2 == 0:
            headers.append(f" {12.0 * factor} {13.0 * factor} {14.0}")
        else:
            headers.append(
                f" {12.0 * factor} {13.0 * factor} {14.0} 80.0 95.0 103.0"
            )

    return _write_trajectory(path, positions, headers)


def _run_rdf(traj, **kwargs):
    """
    Runs an RDF analysis and returns its stacked output arrays.
    """
    analysis = RDF(traj, "O", "H", **kwargs)
    results = np.column_stack(analysis.run())

    return analysis, results


class TestRDFKernelEquivalence:

    """
    Equivalence tests of the Cython kernel vs the numpy fallback.
    """

    def test_orthorhombic(self):
        _assert_kernel_equivalence([Cell(11.0, 13.0, 17.0)] * 40)

    def test_orthorhombic_r_min(self):
        _assert_kernel_equivalence(
            [Cell(11.0, 13.0, 17.0)] * 40,
            r_min=1.5,
            delta_r=0.037,
            n_bins=150,
        )

    def test_triclinic(self):
        _assert_kernel_equivalence(
            [Cell(12.0, 14.0, 16.0, 80.0, 95.0, 103.0)] * 40
        )

    def test_npt_changing_boxes(self):
        # box breathes every 5 frames (NPT-like), including
        # triclinic stretches
        cells = []
        for i in range(40):
            factor = 1.0 + 0.02 * ((i // 5) % 4)
            if (i // 5) % 2 == 0:
                cells.append(
                    Cell(11.0 * factor, 13.0 * factor, 17.0 * factor)
                )
            else:
                cells.append(
                    Cell(
                        11.0 * factor,
                        13.0 * factor,
                        17.0 * factor,
                        85.0,
                        92.0,
                        88.0,
                    )
                )

        _assert_kernel_equivalence(cells)

    def test_vacuum(self):
        _assert_kernel_equivalence([Cell()] * 40, delta_r=0.5, n_bins=60)

    def test_distances_on_bin_edges(self):
        # positions on an exact grid produce distances exactly on
        # the bin edges, exercising the exact floor-divide path
        if _rdf_kernel is None:  # pragma: no cover - build-dependent
            pytest.skip("Cython _rdf_kernel extension not built")

        frames = [
            np.array(
                [[0.25 * i, 0.0, 0.0] for i in range(30)],
                dtype=np.float32,
            )
        ]
        cells = [Cell(20.0, 20.0, 20.0)]

        histograms = [
            _drive_kernel(kernel, frames, cells, 0.0, 0.25, 40)
            for kernel in (
                _rdf_kernel.rdf_frame_histogram,
                _rdf_kernel_py.rdf_frame_histogram,
            )
        ]

        assert np.array_equal(histograms[0], histograms[1])
        assert histograms[0].sum() > 0


class TestRDFFastPath:

    """
    Bit-identity tests of the raw-frame fast path against the
    original RDF path (with both kernel implementations).
    """

    @pytest.mark.parametrize("kernel", KERNELS)
    @pytest.mark.parametrize(
        "write_file", [_write_random_trajectory, _write_npt_trajectory]
    )
    def test_fast_path_matches_original_paths(
        self,
        kernel,
        write_file,
        tmp_path,
        monkeypatch,
    ):
        # the fast path (with either kernel implementation) must
        # reproduce the results of the original path bit for bit,
        # both from an in-memory trajectory and from a reader with
        # the fast path disabled
        filename = write_file(tmp_path / "traj.xyz")

        monkeypatch.setattr(rdf_module, "rdf_frame_histogram", kernel)

        rdf_fast, results_fast = _run_rdf(
            TrajectoryReader(filename), delta_r=0.1, r_max=5.0
        )

        assert rdf_fast._raw_reader is not None

        rdf_memory, results_memory = _run_rdf(
            TrajectoryReader(filename).read(), delta_r=0.1, r_max=5.0
        )

        assert rdf_memory._raw_reader is None

        monkeypatch.setattr(
            RDF, "_use_raw_fast_path", lambda self, traj: False
        )

        rdf_reader, results_reader = _run_rdf(
            TrajectoryReader(filename), delta_r=0.1, r_max=5.0
        )

        assert rdf_reader._raw_reader is None

        assert np.array_equal(rdf_fast.bins, rdf_memory.bins)
        assert np.array_equal(results_fast, results_memory)
        assert np.array_equal(results_fast, results_reader)
        assert rdf_fast.bins.sum() > 0

    @pytest.mark.parametrize("kernel", KERNELS)
    def test_fast_path_matches_original_path_r_min(
        self,
        kernel,
        tmp_path,
        monkeypatch,
    ):
        filename = _write_random_trajectory(tmp_path / "traj.xyz")

        monkeypatch.setattr(rdf_module, "rdf_frame_histogram", kernel)

        rdf_fast, results_fast = _run_rdf(
            TrajectoryReader(filename), delta_r=0.1, r_max=5.0, r_min=1.0
        )

        assert rdf_fast._raw_reader is not None

        rdf_memory, results_memory = _run_rdf(
            TrajectoryReader(filename).read(),
            delta_r=0.1,
            r_max=5.0,
            r_min=1.0,
        )

        assert np.array_equal(rdf_fast.bins, rdf_memory.bins)
        assert np.array_equal(results_fast, results_memory)

    def test_fast_path_multiple_files(self, tmp_path):
        filename1 = _write_random_trajectory(tmp_path / "traj1.xyz", seed=1)
        filename2 = _write_random_trajectory(tmp_path / "traj2.xyz", seed=2)

        rdf_fast, results_fast = _run_rdf(
            TrajectoryReader([filename1, filename2]),
            delta_r=0.1,
            r_max=5.0,
        )

        assert rdf_fast._raw_reader is not None
        assert rdf_fast.n_frames == 50

        traj = TrajectoryReader([filename1, filename2]).read()
        _rdf_memory, results_memory = _run_rdf(
            traj, delta_r=0.1, r_max=5.0
        )

        assert np.array_equal(results_fast, results_memory)

    def test_dispatch_restrictions(self, tmp_path):
        filename = _write_random_trajectory(tmp_path / "traj.xyz")

        # plain reader input uses the fast path
        rdf_fast = RDF(
            TrajectoryReader(filename), "O", "H", delta_r=0.1, r_max=5.0
        )
        assert rdf_fast._raw_reader is not None
        assert rdf_fast.frame_generator is None

        # intra-molecular exclusion keeps the original path
        rdf_no_intra = RDF(
            TrajectoryReader(filename),
            "O",
            "H",
            delta_r=0.1,
            r_max=5.0,
            no_intra_molecular=True,
        )
        assert rdf_no_intra._raw_reader is None
        assert rdf_no_intra.frame_generator is not None

        # in-memory trajectories keep the original path
        rdf_memory = RDF(
            TrajectoryReader(filename).read(),
            "O",
            "H",
            delta_r=0.1,
            r_max=5.0,
        )
        assert rdf_memory._raw_reader is None
        assert rdf_memory.frame_generator is not None

    def test_fast_path_frame_with_missing_atoms(
        self, tmp_path, caplog, monkeypatch
    ):
        filename = _write_random_trajectory(tmp_path / "traj.xyz")

        rdf_fast = RDF(
            TrajectoryReader(filename), "O", "H", delta_r=0.1, r_max=5.0
        )

        def _truncated_generator():
            yield (
                np.zeros((3, 3), dtype=np.float32),
                Cell(12.0, 12.0, 12.0),
            )

        monkeypatch.setattr(
            rdf_fast._raw_reader,
            "raw_frame_generator",
            _truncated_generator,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
                "Frame 1 of the trajectory provides only 3 atoms, "
                "but the selections reference the atom index 15. "
                "Please provide a trajectory with a consistent "
                "number of atoms."
            ),
            exception=RDFError,
            function=rdf_fast.run,
        )

    def test_active_kernel_is_a_known_implementation(self):
        # the rdf module must have wired up either the Cython kernel
        # or the numpy fallback via the try-import
        assert rdf_module.rdf_frame_histogram.__module__ in (
            "PQAnalysis.analysis.rdf._rdf_kernel",
            "PQAnalysis.analysis.rdf._rdf_kernel_py",
        )


class TestScanCells:

    """
    Tests of the header-only cell scan of the fast path against the
    cells full scan of the TrajectoryReader.
    """

    @staticmethod
    def _assert_cells_match(scanned, reference):
        assert len(scanned) == len(reference)

        for scanned_cell, reference_cell in zip(scanned, reference):
            assert np.array_equal(
                scanned_cell.box_lengths, reference_cell.box_lengths
            )
            assert np.array_equal(
                scanned_cell.box_angles, reference_cell.box_angles
            )

    def test_matches_trajectory_reader_cells_npt(self, tmp_path):
        filename = _write_npt_trajectory(tmp_path / "traj.xyz")

        cells, unique_cells = RDF._scan_cells([filename])

        self._assert_cells_match(cells, TrajectoryReader(filename).cells)

        # every distinct box appears exactly once in the unique list
        assert len({id(cell) for cell in cells}) == len(unique_cells)

    def test_matches_trajectory_reader_cells_dedup(self, tmp_path):
        filename = _write_random_trajectory(tmp_path / "traj.xyz")

        cells, unique_cells = RDF._scan_cells([filename])

        self._assert_cells_match(cells, TrajectoryReader(filename).cells)

        assert len(unique_cells) == 1
        assert all(cell is unique_cells[0] for cell in cells)

    def test_matches_trajectory_reader_cells_vacuum_inheritance(
        self, tmp_path
    ):
        rng = np.random.default_rng(11)
        positions = rng.uniform(0.0, 10.0, size=(6, 4, 3))

        # frames without box information inherit the last box
        headers = [
            " 10.0 10.0 10.0",
            "",
            " 11.0 11.0 11.0 80.0 95.0 103.0",
            "",
            "",
            " 10.0 10.0 10.0",
        ]
        filename = _write_trajectory(tmp_path / "traj.xyz", positions, headers)

        cells, unique_cells = RDF._scan_cells([filename])

        self._assert_cells_match(cells, TrajectoryReader(filename).cells)

        assert len(unique_cells) == 2
        assert cells[1] is cells[0]
        assert cells[3] is cells[2]
        assert cells[4] is cells[2]
        assert cells[5] is cells[0]

    def test_matches_trajectory_reader_cells_pure_vacuum(self, tmp_path):
        rng = np.random.default_rng(12)
        positions = rng.uniform(0.0, 10.0, size=(3, 4, 3))

        filename = _write_trajectory(
            tmp_path / "traj.xyz", positions, ["", "", ""]
        )

        cells, unique_cells = RDF._scan_cells([filename])

        self._assert_cells_match(cells, TrajectoryReader(filename).cells)

        assert len(unique_cells) == 1
        assert unique_cells[0].is_vacuum

    def test_multiple_files_inherit_across_boundaries(self, tmp_path):
        rng = np.random.default_rng(13)
        positions = rng.uniform(0.0, 10.0, size=(2, 4, 3))

        filename1 = _write_trajectory(
            tmp_path / "traj1.xyz", positions, [" 10.0 10.0 10.0", ""]
        )
        filename2 = _write_trajectory(
            tmp_path / "traj2.xyz", positions, ["", " 11.0 11.0 11.0"]
        )

        cells, unique_cells = RDF._scan_cells([filename1, filename2])

        reader = TrajectoryReader([filename1, filename2])
        self._assert_cells_match(cells, reader.cells)

        assert len(cells) == 4
        assert cells[2] is cells[0]
        assert len(unique_cells) == 2

    def test_invalid_box_header(self, tmp_path, caplog):
        rng = np.random.default_rng(14)
        positions = rng.uniform(0.0, 10.0, size=(2, 4, 3))

        filename = _write_trajectory(
            tmp_path / "traj.xyz",
            positions,
            [" 10.0 10.0 10.0", " 10.0 10.0"],
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=RDF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
                "Invalid number of arguments for box: 3 encountered "
                f"in file {filename}:2 = 4 10.0 10.0"
            ),
            exception=TrajectoryReaderError,
            function=RDF._scan_cells,
            filenames=[filename],
        )
