"""
Equivalence tests of the Cython VACF kernels against the pure
Python/numpy fallback kernels.

Both kernel implementations expose the identical signatures
(``accumulate_frame``, ``weight_frame`` and ``parse_charge_lines``)
and are driven over synthetic velocity/charge data. The per-origin
dot products of ``accumulate_frame`` are accumulated with different
(but both deterministic) float64 summation orders, so the resulting
correlation accumulators must agree to within tight floating point
tolerances; ``weight_frame`` and ``parse_charge_lines`` must be
bitwise identical. The fast path of the VACF class itself is
additionally run against both kernel implementations via
monkeypatching and compared to the in-memory Trajectory path.
"""

import sys

import numpy as np
import pytest

from PQAnalysis.analysis.vacf import VACF
from PQAnalysis.analysis.vacf import _vacf_kernel_py
from PQAnalysis.analysis.vacf import _raw_charge_reader
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.traj import TrajectoryFormat

from .. import pytestmark  # pylint: disable=unused-import

# pylint: disable=protected-access

try:
    from PQAnalysis.analysis.vacf import _vacf_kernel
except ModuleNotFoundError:  # pragma: no cover - build-dependent
    _vacf_kernel = None

#: The module defining the VACF class (the package attribute ``vacf``
#: is shadowed by the api function of the same name).
vacf_module = sys.modules[VACF.__module__]

KERNEL_MODULES = [
    pytest.param(_vacf_kernel_py, id="python-fallback"),
    pytest.param(
        _vacf_kernel,
        id="cython",
        marks=pytest.mark.skipif(
            _vacf_kernel is None,
            reason="Cython _vacf_kernel extension not built",
        ),
    ),
]



def _require_cython():
    """
    Skips the calling test when the Cython extension is not built.
    """
    if _vacf_kernel is None:  # pragma: no cover - build-dependent
        pytest.skip("Cython _vacf_kernel extension not built")



def _float32_velocities(n_frames, n_atoms, seed):
    """
    Builds float64 per-frame velocities that went through the float32
    parse of the trajectory readers.
    """
    rng = np.random.default_rng(seed)

    values = rng.standard_normal((n_frames, n_atoms, 3))
    values = np.asarray(values, dtype=np.float32)

    return [
        np.ascontiguousarray(frame, dtype=np.float64) for frame in values
    ]



def _drive_accumulate(accumulate, velocities, window, gap):
    """
    Drives an ``accumulate_frame`` implementation frame by frame with
    the exact spawn/shift/drain scheduling of ``VACF._run_direct`` and
    returns the full final state.
    """
    n_frames = len(velocities)
    n_target = velocities[0].shape[0]

    stop_frame = (n_frames - window) // gap * gap

    if stop_frame == 0:
        stop_frame = 1

    n_slots = window // gap + 1

    corr = np.zeros(window + 1, dtype=np.float64)
    origin_vel = np.zeros((n_slots, n_target, 3), dtype=np.float64)
    origin_norm = np.zeros(n_slots, dtype=np.float64)
    origin_frame = np.zeros(n_slots, dtype=np.longlong)
    n_active = 0
    n_origins = 0

    for frame_number, vel in enumerate(velocities, 1):
        spawn = frame_number % gap == 0 and frame_number <= stop_frame

        n_active = accumulate(
            corr,
            origin_vel,
            origin_norm,
            origin_frame,
            n_active,
            vel,
            frame_number,
            spawn,
            window,
        )

        assert n_active >= 0

        if spawn:
            n_origins += 1

    return {
        "corr": corr / n_origins,
        "origin_vel": origin_vel,
        "origin_norm": origin_norm,
        "origin_frame": origin_frame,
        "n_active": n_active,
        "n_origins": n_origins,
    }



def _assert_accumulate_equivalence(n_frames, n_atoms, window, gap, seed=42):
    """
    Runs both ``accumulate_frame`` implementations over the same
    velocity stream and asserts that the complete final state agrees.
    """
    _require_cython()

    velocities = _float32_velocities(n_frames, n_atoms, seed)

    results = [
        _drive_accumulate(module.accumulate_frame, velocities, window, gap)
        for module in (_vacf_kernel, _vacf_kernel_py)
    ]

    for key in ("corr", "origin_norm"):
        assert np.allclose(
            results[0][key],
            results[1][key],
            rtol=1e-14,
            atol=1e-14,
        ), f"kernel/fallback mismatch in {key}"

    # the origin bookkeeping is bitwise identical (plain copies/shifts)
    assert np.array_equal(results[0]["origin_vel"], results[1]["origin_vel"])
    assert np.array_equal(
        results[0]["origin_frame"],
        results[1]["origin_frame"],
    )
    assert results[0]["n_active"] == results[1]["n_active"]
    assert results[0]["n_origins"] == results[1]["n_origins"]

    # both kernels must have accumulated a non-trivial correlation
    assert results[0]["corr"][0] == 1.0
    assert np.any(results[0]["corr"][1:] != 0.0)



class TestAccumulateFrameEquivalence:

    """
    Equivalence tests of the Cython ``accumulate_frame`` kernel vs the
    numpy fallback.
    """

    def test_gap_one(self):
        _assert_accumulate_equivalence(60, 5, window=12, gap=1)

    def test_gap_larger_one(self):
        _assert_accumulate_equivalence(83, 7, window=20, gap=5)

    def test_single_atom(self):
        _assert_accumulate_equivalence(50, 1, window=10, gap=2)

    def test_window_equals_n_frames_boundary(self):
        # legacy stop_frame reset 0 -> 1: a single origin at frame 1
        # and the final lag bin stays zero
        _assert_accumulate_equivalence(20, 3, window=20, gap=1)

    def test_drain_phase(self):
        # trailing frames past the last origin only drain old origins
        _assert_accumulate_equivalence(47, 4, window=8, gap=4)

    @pytest.mark.parametrize("kernel_module", KERNEL_MODULES)
    def test_zero_norm_spawn_leaves_state_unchanged(self, kernel_module):
        # a zero-norm origin returns -1 and must not modify any state
        corr = np.zeros(6, dtype=np.float64)
        origin_vel = np.zeros((6, 2, 3), dtype=np.float64)
        origin_norm = np.zeros(6, dtype=np.float64)
        origin_frame = np.zeros(6, dtype=np.longlong)

        result = kernel_module.accumulate_frame(
            corr,
            origin_vel,
            origin_norm,
            origin_frame,
            0,
            np.zeros((2, 3), dtype=np.float64),
            1,
            True,
            5,
        )

        assert result == -1
        assert not np.any(corr)
        assert not np.any(origin_vel)
        assert not np.any(origin_norm)
        assert not np.any(origin_frame)



class TestWeightFrameEquivalence:

    """
    Bitwise equivalence tests of the Cython ``weight_frame`` kernel vs
    the numpy fallback.
    """

    @staticmethod
    def _raw_values(n_atoms=9, seed=7):
        rng = np.random.default_rng(seed)

        return np.asarray(
            rng.standard_normal((n_atoms, 3)),
            dtype=np.float32,
        )

    def test_without_charges(self):
        _require_cython()

        values = self._raw_values()
        indices = np.array([0, 2, 3, 8], dtype=np.intp)

        result = _vacf_kernel.weight_frame(values, indices, None)
        reference = _vacf_kernel_py.weight_frame(values, indices, None)

        assert result.dtype == np.float64
        assert np.array_equal(result, reference)
        # the float32 -> float64 cast is exact
        assert np.array_equal(
            result,
            np.asarray(values, dtype=np.float64)[indices],
        )

    def test_with_charges(self):
        _require_cython()

        values = self._raw_values(seed=11)
        indices = np.array([1, 4, 5, 6, 7], dtype=np.intp)
        charges = np.array([-0.8, 0.4, 0.4, -1.2, 2.0], dtype=np.float64)

        result = _vacf_kernel.weight_frame(values, indices, charges)
        reference = _vacf_kernel_py.weight_frame(values, indices, charges)

        assert result.dtype == np.float64
        assert np.array_equal(result, reference)

    def test_qmcfc_stripped_view(self):
        # the QMCFC fast path passes the values without the leading
        # dummy atom row (a row-sliced view of the parsed array)
        _require_cython()

        values = self._raw_values(n_atoms=5, seed=13)[1:]
        indices = np.arange(4, dtype=np.intp)
        charges = np.array([0.1, -0.2, 0.3, -0.4], dtype=np.float64)

        assert np.array_equal(
            _vacf_kernel.weight_frame(values, indices, charges),
            _vacf_kernel_py.weight_frame(values, indices, charges),
        )



class TestParseChargeLinesEquivalence:

    """
    Equivalence tests of the Cython ``parse_charge_lines`` kernel vs
    the numpy fallback.
    """

    VALID_LINES = [
        "O -0.89076318\n",
        "H 0.44538159\n",
        "  Na1 1.0  \n",
        "Cl\t-1.25e-03\n",
        "X 0.0\n",
        "O2 .5\n",
        "H +4.75\n",
        "C -17\n",
    ]

    INVALID_LINES = [
        "O\n",
        "O 1.0 2.0\n",
        "O abc\n",
        "O 1.0x\n",
        "O 1,5\n",
        "\n",
    ]

    def test_valid_lines_bitwise_identical(self):
        _require_cython()

        n_atoms = len(self.VALID_LINES)

        result = _vacf_kernel.parse_charge_lines(self.VALID_LINES, n_atoms)
        reference = _vacf_kernel_py.parse_charge_lines(
            self.VALID_LINES,
            n_atoms,
        )

        assert result.dtype == np.float64
        assert np.array_equal(result, reference)
        assert np.array_equal(
            result,
            np.array(
                [float(line.split()[1]) for line in self.VALID_LINES]
            ),
        )

    @pytest.mark.parametrize("kernel_module", KERNEL_MODULES)
    @pytest.mark.parametrize(
        "bad_line",
        INVALID_LINES,
        ids=repr,
    )
    def test_invalid_lines_raise(self, kernel_module, bad_line):
        with pytest.raises(ValueError):
            kernel_module.parse_charge_lines(["O 1.0\n", bad_line], 2)



class TestVACFFastPathKernels:

    """
    Tests of the VACF fast path with both kernel implementations.
    """

    N_FRAMES = 60
    N_ATOMS = 6

    @classmethod
    def _write_trajectories(cls, tmp_path, seed=2026):
        """
        Writes a velocity and a lockstep charge trajectory file and
        returns their paths.
        """
        rng = np.random.default_rng(seed)

        velocities = rng.standard_normal((cls.N_FRAMES, cls.N_ATOMS, 3))
        charges = 0.5 + 0.1 * rng.standard_normal(
            (cls.N_FRAMES, cls.N_ATOMS)
        )

        names = [
            "O" if i % 2 == 0 else "H" for i in range(cls.N_ATOMS)
        ]

        vel_file = str(tmp_path / "traj.vel")
        charge_file = str(tmp_path / "traj.chrg")

        with open(vel_file, "w", encoding="utf-8") as file:
            for frame in velocities:
                file.write(f"{cls.N_ATOMS} 10.0 11.0 12.0\n\n")
                for name, (x, y, z) in zip(names, frame):
                    file.write(f"{name} {x:.7f} {y:.7f} {z:.7f}\n")

        with open(charge_file, "w", encoding="utf-8") as file:
            for frame in charges:
                file.write(f"{cls.N_ATOMS} 10.0 11.0 12.0\n\n")
                for name, charge in zip(names, frame):
                    file.write(f"{name} {charge:.7f}\n")

        return vel_file, charge_file

    def _patch_kernels(self, monkeypatch, kernel_module):
        monkeypatch.setattr(
            vacf_module,
            "accumulate_frame",
            kernel_module.accumulate_frame,
        )
        monkeypatch.setattr(
            vacf_module,
            "weight_frame",
            kernel_module.weight_frame,
        )
        monkeypatch.setattr(
            _raw_charge_reader,
            "parse_charge_lines",
            kernel_module.parse_charge_lines,
        )

    @staticmethod
    def _charge_kwargs(charge_source, vel_file, charge_file):
        if charge_source == "static":
            return {"charges": np.array([-0.8, 0.4] * 3)}

        if charge_source == "trajectory":
            return {
                "charge_traj": TrajectoryReader(
                    charge_file,
                    traj_format=TrajectoryFormat.CHARGE,
                )
            }

        return {}

    @pytest.mark.parametrize(
        "charge_source",
        ["none", "static", "trajectory"],
    )
    @pytest.mark.parametrize("kernel_module", KERNEL_MODULES)
    def test_fast_path_matches_in_memory_path(
        self,
        kernel_module,
        charge_source,
        tmp_path,
        monkeypatch,
    ):
        # the fast path (with either kernel implementation) must
        # reproduce the results of the in-memory Trajectory hot loop
        # for all charge modes
        vel_file, charge_file = self._write_trajectories(tmp_path)

        self._patch_kernels(monkeypatch, kernel_module)

        fast = VACF(
            TrajectoryReader(vel_file),
            window_size=20,
            time_step=0.1,
            gap=5,
            **self._charge_kwargs(charge_source, vel_file, charge_file),
        )

        assert fast._raw_reader is not None

        if charge_source == "trajectory":
            assert fast._raw_charge_reader is not None

        _, fast_correlation = fast.run()

        charge_kwargs = self._charge_kwargs(
            charge_source,
            vel_file,
            charge_file,
        )

        if charge_source == "trajectory":
            charge_kwargs = {
                "charge_traj":
                TrajectoryReader(
                    charge_file,
                    traj_format=TrajectoryFormat.CHARGE,
                ).read()
            }

        reference = VACF(
            TrajectoryReader(vel_file).read(),
            window_size=20,
            time_step=0.1,
            gap=5,
            **charge_kwargs,
        )

        assert reference._raw_reader is None

        _, reference_correlation = reference.run()

        assert fast.n_origins == reference.n_origins
        assert np.allclose(
            fast_correlation,
            reference_correlation,
            rtol=0.0,
            atol=1e-14,
        )

    @pytest.mark.parametrize("kernel_module", KERNEL_MODULES)
    def test_fast_path_matches_in_memory_path_fft(
        self,
        kernel_module,
        tmp_path,
        monkeypatch,
    ):
        # the fft estimator consumes the same velocity stream
        vel_file, _ = self._write_trajectories(tmp_path)

        self._patch_kernels(monkeypatch, kernel_module)

        fast = VACF(
            TrajectoryReader(vel_file),
            window_size=20,
            time_step=0.1,
            method="fft",
        )

        assert fast._raw_reader is not None

        _, fast_correlation = fast.run()

        reference = VACF(
            TrajectoryReader(vel_file).read(),
            window_size=20,
            time_step=0.1,
            method="fft",
        )

        _, reference_correlation = reference.run()

        assert np.allclose(
            fast_correlation,
            reference_correlation,
            rtol=0.0,
            atol=1e-14,
        )

    def test_active_kernels_are_a_known_implementation(self):
        # the vacf modules must have wired up either the Cython kernel
        # or the numpy fallback via the try-import
        known_modules = (
            "PQAnalysis.analysis.vacf._vacf_kernel",
            "PQAnalysis.analysis.vacf._vacf_kernel_py",
        )

        assert vacf_module.accumulate_frame.__module__ in known_modules
        assert vacf_module.weight_frame.__module__ in known_modules
        assert (
            _raw_charge_reader.parse_charge_lines.__module__
            in known_modules
        )
