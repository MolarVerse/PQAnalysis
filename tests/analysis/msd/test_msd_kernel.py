"""
Equivalence tests of the Cython MSD frame kernel against the pure
Python/numpy fallback kernel.

Both kernels implement the identical signature and are driven frame by
frame over synthetic wrapped random-walk trajectories (orthorhombic,
triclinic, NPT with changing boxes and vacuum cases). The resulting
MSD accumulators and the complete mutable state (shift, prev_pos,
origins, bookkeeping) must agree to within tight floating point
tolerances. The fast path of the MSD class itself is additionally run
against both kernels via monkeypatching.
"""

import sys

import numpy as np
import pytest

from PQAnalysis.analysis.msd import MSD
from PQAnalysis.analysis.msd import _msd_kernel_py
from PQAnalysis.core import Cell
from PQAnalysis.io import TrajectoryReader

from .. import pytestmark  # pylint: disable=unused-import

try:
    from PQAnalysis.analysis.msd import _msd_kernel
except ModuleNotFoundError:  # pragma: no cover - build-dependent
    _msd_kernel = None

#: The module defining the MSD class (the package attribute ``msd``
#: is shadowed by the api function of the same name).
msd_module = sys.modules[MSD.__module__]

KERNELS = [
    pytest.param(_msd_kernel_py.msd_frame_update, id="python-fallback"),
    pytest.param(
        _msd_kernel.msd_frame_update if _msd_kernel is not None else None,
        id="cython",
        marks=pytest.mark.skipif(
            _msd_kernel is None,
            reason="Cython _msd_kernel extension not built",
        ),
    ),
]



def _wrapped_random_walk(n_frames, n_atoms, cells, seed):
    """
    Builds a float32 wrapped random walk (one array per frame) that is
    wrapped into the (possibly per-frame changing) cells.
    """
    rng = np.random.default_rng(seed)

    positions = np.cumsum(
        rng.normal(0.0, 1.5, size=(n_frames, n_atoms, 3)),
        axis=0,
    )

    frames = []

    for frame_positions, cell in zip(positions, cells):
        if cell.is_vacuum:
            wrapped = frame_positions
        else:
            # wrap into the cell with the generic triclinic recipe
            fractional = frame_positions @ np.linalg.inv(cell.box_matrix).T
            wrapped = (
                frame_positions
                - np.floor(fractional) @ cell.box_matrix.T
            )

        frames.append(np.asarray(wrapped, dtype=np.float32))

    return frames



def _drive_kernel(
    kernel,
    frames,
    cells,
    indices,
    window,
    gap,
    n_start,
    stop_frame,
    n_origins_max,
):
    """
    Drives a kernel implementation frame by frame and returns the
    full final state.
    """
    n_sel = len(indices)

    msd = np.zeros((window + 1, 3))
    origins = np.zeros((n_origins_max, n_sel, 3))
    state = np.zeros(2, dtype=np.int64)

    pos = np.zeros((n_sel, 3))
    prev_pos = np.zeros((n_sel, 3))
    shift = np.zeros((n_sel, 3))
    unwrapped = np.zeros((n_sel, 3))

    indices = np.ascontiguousarray(indices, dtype=np.int64)

    for counter, (values, cell) in enumerate(zip(frames, cells), start=1):
        is_vacuum = 1 if cell.is_vacuum else 0

        if is_vacuum:
            box = np.eye(3)
            inv_box = np.eye(3)
        else:
            box = np.ascontiguousarray(cell.box_matrix, dtype=np.float64)
            inv_box = np.ascontiguousarray(
                cell.inverse_box_matrix, dtype=np.float64
            )

        kernel(
            values,
            indices,
            box,
            inv_box,
            is_vacuum,
            pos,
            prev_pos,
            shift,
            unwrapped,
            origins,
            msd,
            state,
            counter,
            gap,
            window,
            n_start,
            stop_frame,
        )

    return {
        "msd": msd,
        "origins": origins,
        "state": state,
        "pos": pos,
        "prev_pos": prev_pos,
        "shift": shift,
        "unwrapped": unwrapped,
    }



def _assert_kernel_equivalence(cells, n_frames, window, gap, n_start=0):
    """
    Runs both kernels over the same trajectory and asserts that the
    complete final state agrees.
    """
    if _msd_kernel is None:  # pragma: no cover - build-dependent
        pytest.skip("Cython _msd_kernel extension not built")

    n_atoms = 7
    frames = _wrapped_random_walk(n_frames, n_atoms, cells, seed=2024)
    indices = np.array([0, 2, 3, 6], dtype=np.int64)

    stop_frame = (n_frames - window) // gap * gap
    n_origins_max = window // gap

    results = [
        _drive_kernel(
            kernel,
            frames,
            cells,
            indices,
            window,
            gap,
            n_start,
            stop_frame,
            n_origins_max,
        )
        for kernel in (
            _msd_kernel.msd_frame_update,
            _msd_kernel_py.msd_frame_update,
        )
    ]

    for key in results[0]:
        assert np.allclose(
            results[0][key],
            results[1][key],
            rtol=0.0,
            atol=1e-12,
        ), f"kernel/fallback mismatch in {key}"

    assert np.array_equal(results[0]["state"], results[1]["state"])

    # both kernels must have accumulated a non-trivial MSD
    assert np.any(results[0]["msd"] != 0.0)



class TestMSDKernelEquivalence:

    """
    Equivalence tests of the Cython kernel vs the numpy fallback.
    """

    def test_orthorhombic(self):
        cells = [Cell(11.0, 13.0, 17.0)] * 90

        _assert_kernel_equivalence(cells, 90, window=20, gap=5)

    def test_orthorhombic_gap_one(self):
        cells = [Cell(9.0, 9.0, 9.0)] * 50

        _assert_kernel_equivalence(cells, 50, window=10, gap=1)

    def test_orthorhombic_n_start(self):
        cells = [Cell(11.0, 13.0, 17.0)] * 90

        _assert_kernel_equivalence(cells, 90, window=20, gap=5, n_start=17)

    def test_triclinic(self):
        cells = [Cell(12.0, 14.0, 16.0, 80.0, 95.0, 103.0)] * 90

        _assert_kernel_equivalence(cells, 90, window=20, gap=5)

    def test_npt_changing_boxes(self):
        # box breathes every 10 frames (NPT-like), including a
        # triclinic stretch
        cells = []
        for i in range(90):
            factor = 1.0 + 0.02 * ((i // 10) % 4)
            if (i // 10) % 2 == 0:
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

        _assert_kernel_equivalence(cells, 90, window=20, gap=5)

    def test_vacuum(self):
        cells = [Cell()] * 60

        _assert_kernel_equivalence(cells, 60, window=20, gap=5)

    def test_mixed_vacuum_and_box(self):
        # vacuum frames in between (no unwrapping applied there)
        cells = [
            Cell() if i % 7 == 3 else Cell(11.0, 13.0, 17.0)
            for i in range(60)
        ]

        _assert_kernel_equivalence(cells, 60, window=20, gap=5)



class TestMSDFastPathKernels:

    """
    Tests of the MSD fast path with both kernel implementations.
    """

    @staticmethod
    def _write_trajectory(path, n_frames=60, n_atoms=8, seed=99):
        rng = np.random.default_rng(seed)
        box = np.array([10.0, 12.0, 14.0])

        positions = np.cumsum(
            rng.normal(0.0, 0.8, size=(n_frames, n_atoms, 3)),
            axis=0,
        ) % box

        names = ["O" if i % 2 == 0 else "H" for i in range(n_atoms)]

        with open(path, "w", encoding="utf-8") as file:
            for frame_positions in positions:
                file.write(f"{n_atoms} {box[0]} {box[1]} {box[2]}\n\n")
                for name, (x, y, z) in zip(names, frame_positions):
                    file.write(f"{name} {x:.6f} {y:.6f} {z:.6f}\n")

        return positions

    @pytest.mark.parametrize("kernel", KERNELS)
    def test_fast_path_matches_in_memory_path(
        self,
        kernel,
        tmp_path,
        monkeypatch,
    ):
        # the fast path (with either kernel implementation) must
        # reproduce the results of the original in-memory hot loop
        filename = str(tmp_path / "traj.xyz")
        self._write_trajectory(filename)

        monkeypatch.setattr(msd_module, "msd_frame_update", kernel)

        reader = TrajectoryReader(filename)
        msd_fast = MSD(reader, "O", window=20, gap=5)

        assert msd_fast._raw_reader is not None  # pylint: disable=protected-access

        result_fast = np.column_stack(msd_fast.run()[1:])

        traj = TrajectoryReader(filename).read()
        msd_reference = MSD(traj, "O", window=20, gap=5)

        assert msd_reference._raw_reader is None  # pylint: disable=protected-access

        result_reference = np.column_stack(msd_reference.run()[1:])

        assert np.allclose(
            result_fast, result_reference, rtol=0.0, atol=1e-12
        )

    def test_active_kernel_is_a_known_implementation(self):
        # the msd module must have wired up either the Cython kernel
        # or the numpy fallback via the try-import
        assert msd_module.msd_frame_update.__module__ in (
            "PQAnalysis.analysis.msd._msd_kernel",
            "PQAnalysis.analysis.msd._msd_kernel_py",
        )
