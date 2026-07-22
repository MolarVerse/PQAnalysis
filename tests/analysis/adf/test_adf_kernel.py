"""
Equivalence tests of the Cython ADF angle-histogram kernel against the
pure Python/numpy fallback kernel and bit-identity tests of the
raw-frame fast path against the original ADF path.

Both kernels implement the identical signature and are driven frame by
frame over synthetic float32 frames (orthorhombic, triclinic, NPT with
changing boxes, vacuum, gated and bin-edge cases). The resulting int64
histograms must agree exactly (atol 0). The fast path of the ADF class
itself is additionally run against both kernels via monkeypatching and
compared bit-identically against the original path.
"""

import sys

import numpy as np
import pytest

from PQAnalysis.analysis import ADF
from PQAnalysis.analysis.adf import _adf_kernel_py
from PQAnalysis.core import Cell
from PQAnalysis.io import TrajectoryReader

from .. import pytestmark  # pylint: disable=unused-import

# pylint: disable=protected-access

try:
    from PQAnalysis.analysis.adf import _adf_kernel
except ModuleNotFoundError:  # pragma: no cover - build-dependent
    _adf_kernel = None

#: The module defining the ADF class (the package attribute ``adf`` is
#: shadowed by the api function of the same name).
adf_module = sys.modules[ADF.__module__]

KERNELS = [
    pytest.param(_adf_kernel_py.adf_frame_histogram, id="python-fallback"),
    pytest.param(
        _adf_kernel.adf_frame_histogram if _adf_kernel is not None else None,
        id="cython",
        marks=pytest.mark.skipif(
            _adf_kernel is None,
            reason="Cython _adf_kernel extension not built",
        ),
    ),
]

NO_GATE = (0, 0.0, np.inf, 0, 0.0, np.inf)


def _random_frames(n_frames, n_atoms, spread, seed):
    """Builds float32 frames of uniformly distributed positions."""
    rng = np.random.default_rng(seed)

    return [
        np.asarray(
            rng.uniform(0.0, spread, size=(n_atoms, 3)),
            dtype=np.float32,
        )
        for _ in range(n_frames)
    ]


def _drive_kernel(kernel, frames, cells, gates, delta_angle=1.0, n_bins=180):
    """Drives a kernel implementation frame by frame."""
    reference_indices = np.arange(0, 20, 2, dtype=np.int64)
    # overlapping index sets exercise the self-pair exclusions
    target_indices = np.arange(0, 20, 3, dtype=np.int64)
    target_indices_2 = np.arange(1, 20, 2, dtype=np.int64)

    hist = np.zeros(n_bins, dtype=np.int64)

    for values, cell in zip(frames, cells):
        is_orthorhombic = 1 if (
            cell.alpha == 90 and cell.beta == 90 and cell.gamma == 90
        ) else 0

        kernel(
            values,
            reference_indices,
            target_indices,
            target_indices_2,
            np.ascontiguousarray(cell.box_lengths, dtype=np.float64),
            np.ascontiguousarray(cell.box_matrix, dtype=np.float64),
            np.ascontiguousarray(cell.inverse_box_matrix, dtype=np.float64),
            is_orthorhombic,
            gates[0],
            gates[1],
            gates[2],
            gates[3],
            gates[4],
            gates[5],
            delta_angle,
            n_bins,
            hist,
        )

    return hist


def _assert_kernel_equivalence(
    cells, gates=NO_GATE, delta_angle=1.0, n_bins=180
):
    """Runs both kernels over the same frames and asserts exact agreement."""
    if _adf_kernel is None:  # pragma: no cover - build-dependent
        pytest.skip("Cython _adf_kernel extension not built")

    frames = _random_frames(len(cells), 20, 15.0, seed=2026)

    histograms = [
        _drive_kernel(kernel, frames, cells, gates, delta_angle, n_bins)
        for kernel in (
            _adf_kernel.adf_frame_histogram,
            _adf_kernel_py.adf_frame_histogram,
        )
    ]

    assert np.array_equal(histograms[0], histograms[1])
    assert histograms[0].sum() > 0



class TestADFKernelEquivalence:

    def test_orthorhombic(self):
        _assert_kernel_equivalence([Cell(15.0, 15.0, 15.0)] * 20)

    def test_orthorhombic_gated(self):
        _assert_kernel_equivalence(
            [Cell(15.0, 15.0, 15.0)] * 20,
            gates=(1, 1.0, 6.0, 1, 1.0, 6.0),
        )

    def test_orthorhombic_asymmetric_gates(self):
        _assert_kernel_equivalence(
            [Cell(13.0, 15.0, 17.0)] * 20,
            gates=(1, 2.0, 7.0, 0, 0.0, np.inf),
        )

    def test_triclinic(self):
        _assert_kernel_equivalence(
            [Cell(15.0, 16.0, 17.0, 80.0, 95.0, 103.0)] * 20
        )

    def test_triclinic_gated(self):
        _assert_kernel_equivalence(
            [Cell(15.0, 16.0, 17.0, 80.0, 95.0, 103.0)] * 20,
            gates=(1, 1.5, 6.5, 1, 2.0, 7.0),
            delta_angle=2.0,
            n_bins=90,
        )

    def test_npt_changing_boxes(self):
        cells = []
        for i in range(20):
            factor = 1.0 + 0.02 * ((i // 5) % 4)
            if (i // 5) % 2 == 0:
                cells.append(
                    Cell(13.0 * factor, 15.0 * factor, 17.0 * factor)
                )
            else:
                cells.append(
                    Cell(
                        13.0 * factor,
                        15.0 * factor,
                        17.0 * factor,
                        85.0,
                        92.0,
                        88.0,
                    )
                )

        _assert_kernel_equivalence(cells)

    def test_vacuum(self):
        _assert_kernel_equivalence([Cell()] * 20)

    def test_angles_on_bin_edges(self):
        # a cubic grid produces many collinear / right-angle triplets
        # whose angles land exactly on bin edges, exercising the floor
        if _adf_kernel is None:  # pragma: no cover - build-dependent
            pytest.skip("Cython _adf_kernel extension not built")

        grid = np.array(
            [[x, y, 0.0] for x in range(5) for y in range(5)],
            dtype=np.float32,
        )
        frames = [grid]
        cells = [Cell(20.0, 20.0, 20.0)]

        histograms = [
            _drive_kernel(kernel, frames, cells, NO_GATE)
            for kernel in (
                _adf_kernel.adf_frame_histogram,
                _adf_kernel_py.adf_frame_histogram,
            )
        ]

        assert np.array_equal(histograms[0], histograms[1])
        assert histograms[0].sum() > 0



class TestADFFastPath:

    @staticmethod
    def _write_random_trajectory(path, n_frames=6, n_atoms=15, seed=42):
        rng = np.random.default_rng(seed)
        with open(path, "w", encoding="utf-8") as file:
            for _ in range(n_frames):
                pos = rng.uniform(0.0, 12.0, size=(n_atoms, 3))
                file.write(f"{n_atoms} 12.0 12.0 12.0\n\n")
                for i, (x, y, z) in enumerate(pos):
                    name = "O" if i % 3 == 0 else "H"
                    file.write(f"{name} {x:.6f} {y:.6f} {z:.6f}\n")
        return str(path)

    @pytest.mark.parametrize("kernel", KERNELS)
    def test_fast_path_matches_original_paths(
        self, kernel, tmp_path, monkeypatch
    ):
        filename = self._write_random_trajectory(tmp_path / "traj.xyz")

        monkeypatch.setattr(adf_module, "adf_frame_histogram", kernel)

        adf_fast = ADF(TrajectoryReader(filename), "O", "H", n_angle_bins=180)
        results_fast = np.column_stack(adf_fast.run())
        assert adf_fast._raw_reader is not None

        adf_memory = ADF(
            TrajectoryReader(filename).read(), "O", "H", n_angle_bins=180
        )
        results_memory = np.column_stack(adf_memory.run())
        assert adf_memory._raw_reader is None

        monkeypatch.setattr(
            ADF, "_use_raw_fast_path", lambda self, traj: False
        )
        adf_reader = ADF(TrajectoryReader(filename), "O", "H", n_angle_bins=180)
        results_reader = np.column_stack(adf_reader.run())
        assert adf_reader._raw_reader is None

        assert np.array_equal(adf_fast.angle_bins, adf_memory.angle_bins)
        assert np.array_equal(results_fast, results_memory)
        assert np.array_equal(results_fast, results_reader)
        assert adf_fast.angle_bins.sum() > 0

    @pytest.mark.parametrize("kernel", KERNELS)
    def test_fast_path_matches_original_paths_gated(
        self, kernel, tmp_path, monkeypatch
    ):
        filename = self._write_random_trajectory(tmp_path / "traj.xyz")

        monkeypatch.setattr(adf_module, "adf_frame_histogram", kernel)

        kwargs = dict(
            delta_angle=2.0,
            r_min_1=0.8,
            r_max_1=4.0,
            r_min_2=0.8,
            r_max_2=4.0,
        )

        adf_fast = ADF(TrajectoryReader(filename), "O", "H", **kwargs)
        results_fast = np.column_stack(adf_fast.run())

        adf_memory = ADF(
            TrajectoryReader(filename).read(), "O", "H", **kwargs
        )
        results_memory = np.column_stack(adf_memory.run())

        assert np.array_equal(results_fast, results_memory)
        assert adf_fast.angle_bins.sum() > 0

    def test_fast_path_multiple_files(self, tmp_path):
        filename1 = self._write_random_trajectory(
            tmp_path / "traj1.xyz", seed=1
        )
        filename2 = self._write_random_trajectory(
            tmp_path / "traj2.xyz", seed=2
        )

        adf_fast = ADF(
            TrajectoryReader([filename1, filename2]),
            "O",
            "H",
            n_angle_bins=180,
        )
        results_fast = np.column_stack(adf_fast.run())

        assert adf_fast._raw_reader is not None
        assert adf_fast.n_frames == 12

        traj = TrajectoryReader([filename1, filename2]).read()
        adf_memory = ADF(traj, "O", "H", n_angle_bins=180)
        results_memory = np.column_stack(adf_memory.run())

        assert np.array_equal(results_fast, results_memory)

    def test_active_kernel_is_a_known_implementation(self):
        assert adf_module.adf_frame_histogram.__module__ in (
            "PQAnalysis.analysis.adf._adf_kernel",
            "PQAnalysis.analysis.adf._adf_kernel_py",
        )
