"""Full-precision parity tests against the corrected legacy RDF tool.

The oracle hashes were generated from ``thh_tools/RDF`` compiled with
``-O2 -ffp-contract=off``. Three undefined behaviors were corrected in
the oracle only: the initial histogram allocation size, the ``%lf``
header destinations and the uninitialized integration accumulator.
The distance, binning and normalization expressions are unchanged.
"""

import hashlib
import importlib
import math
from pathlib import Path

import numpy as np
import pytest

from PQAnalysis.analysis.rdf import RDF, rdf
from PQAnalysis.analysis.rdf._rdf_kernel_py import (
    legacy_rdf_batch_histogram as legacy_rdf_batch_histogram_py,
    legacy_rdf_frame_histogram as legacy_rdf_frame_histogram_py,
)
from PQAnalysis.io import TrajectoryReader

from .. import pytestmark  # pylint: disable=unused-import

rdf_module = importlib.import_module("PQAnalysis.analysis.rdf.rdf")

try:
    from PQAnalysis.analysis.rdf._rdf_kernel import (  # pylint: disable=import-error
        legacy_rdf_batch_histogram as legacy_rdf_batch_histogram_compiled,
        legacy_rdf_frame_histogram as legacy_rdf_frame_histogram_compiled,
    )
except ModuleNotFoundError:
    legacy_rdf_batch_histogram_compiled = None
    legacy_rdf_frame_histogram_compiled = None


RDF_FLOAT64_SHA256 = (
    "025ff1e0d08148e5bf5d358ea3318143ff7dc047daa62bd5c2e2119b2e4cf8e4"
)
RDF_HISTOGRAM_SHA256 = (
    "80acbcf723e83519c568f63be37d712a5ff0ae6afce9352f900d4aaf949c3499"
)

KERNELS = [(legacy_rdf_frame_histogram_py, legacy_rdf_batch_histogram_py)]

if legacy_rdf_frame_histogram_compiled is not None:
    KERNELS.append(
        (
            legacy_rdf_frame_histogram_compiled,
            legacy_rdf_batch_histogram_compiled,
        )
    )


def _float64_sha256(values):
    array = np.ascontiguousarray(values, dtype="<f8")

    return hashlib.sha256(array.tobytes()).hexdigest()


def _int64_sha256(values):
    array = np.ascontiguousarray(values, dtype="<i8")

    return hashlib.sha256(array.tobytes()).hexdigest()


@pytest.mark.parametrize("example_dir", ["msd"], indirect=False)
@pytest.mark.parametrize(("frame_kernel", "batch_kernel"), KERNELS)
def test_full_precision_legacy_parity(
    test_with_data_dir,
    monkeypatch,
    frame_kernel,
    batch_kernel,
):
    monkeypatch.setattr(
        rdf_module, "legacy_rdf_frame_histogram", frame_kernel
    )
    monkeypatch.setattr(
        rdf_module, "legacy_rdf_batch_histogram", batch_kernel
    )

    analysis = RDF(
        TrajectoryReader("traj.xyz"),
        "O",
        "H",
        delta_r=0.1,
    )
    result = analysis.run()

    assert analysis._legacy_rdf  # pylint: disable=protected-access
    assert analysis._raw_reader.dtype == "float64"  # pylint: disable=protected-access
    assert analysis.n_bins == 99
    assert analysis.delta_r == float(np.float32(0.1))
    assert _int64_sha256(analysis.bins) == RDF_HISTOGRAM_SHA256
    assert _float64_sha256(np.column_stack(result)) == RDF_FLOAT64_SHA256


@pytest.mark.parametrize(("frame_kernel", "batch_kernel"), KERNELS)
def test_float64_coordinate_controls_boundary_bin(
    tmp_path,
    monkeypatch,
    frame_kernel,
    batch_kernel,
):
    monkeypatch.setattr(
        rdf_module, "legacy_rdf_frame_histogram", frame_kernel
    )
    monkeypatch.setattr(
        rdf_module, "legacy_rdf_batch_histogram", batch_kernel
    )

    trajectory = tmp_path / "boundary.xyz"
    trajectory.write_text(
        "2 1.0 1.0 1.0\n\n"
        "X 0.0 0.0 0.0\n"
        "Y 0.124999999 0.0 0.0\n",
        encoding="utf-8",
    )

    analysis = RDF(
        TrajectoryReader(str(trajectory)),
        "X",
        "Y",
        delta_r=0.125,
    )
    analysis.run()

    assert np.array_equal(analysis.bins, np.array([1.0, 0.0, 0.0, 0.0]))


@pytest.mark.parametrize(("frame_kernel", "batch_kernel"), KERNELS)
def test_noncubic_box_bins_span_only_the_kernel_cutoff(
    tmp_path,
    monkeypatch,
    frame_kernel,
    batch_kernel,
):
    monkeypatch.setattr(
        rdf_module, "legacy_rdf_frame_histogram", frame_kernel
    )
    monkeypatch.setattr(
        rdf_module, "legacy_rdf_batch_histogram", batch_kernel
    )

    trajectory = tmp_path / "noncubic.xyz"
    trajectory.write_text(
        "2 20.0 20.0 60.0 90.0 90.0 90.0\n\n"
        "X 0.0 0.0 0.0\n"
        "Y 3.0 0.0 0.0\n",
        encoding="utf-8",
    )

    analysis = RDF(
        TrajectoryReader(str(trajectory)),
        "X",
        "Y",
        delta_r=0.1,
    )
    result = analysis.run()

    # The kernel drops every pair beyond min(20, 20, 60) / 2 = 10.0, so
    # the histogram must cover int(10.0 / float32(0.1)) = 99 bins. Sizing
    # it from the largest edge (299 bins) would normalize over a tail the
    # kernel never counts.
    delta_r = float(np.float32(0.1))

    assert analysis._legacy_rdf  # pylint: disable=protected-access
    assert analysis.n_bins == 99
    assert analysis.r_max == 99 * delta_r

    # The single pair sits at distance 3.0 -> bin floor(3.0 / delta_r) = 29.
    expected_bins = np.zeros(99)
    expected_bins[29] = 1.0
    assert np.array_equal(analysis.bins, expected_bins)

    # Hand-computed: one pair, one reference, one frame, number density
    # 1 / (20 * 20 * 60), ideal-gas shell between 29 and 30 bin edges.
    density = 1.0 / (20.0 * 20.0 * 60.0)
    shell = 4.0 / 3.0 * math.pi * (
        (30.0 * delta_r) ** 3 - (29.0 * delta_r) ** 3
    )
    assert result[1][29] == pytest.approx(1.0 / (density * shell), rel=1e-12)
    assert result[2][-1] == 1.0


@pytest.mark.parametrize("example_dir", ["msd"], indirect=False)
def test_api_output_round_trips_exact_values(test_with_data_dir):
    Path("rdf.in").write_text(
        "traj_files = traj.xyz\n"
        "reference_selection = O\n"
        "target_selection = H\n"
        "out_file = rdf.dat\n"
        "delta_r = 0.1\n",
        encoding="utf-8",
    )

    rdf("rdf.in")

    result = np.loadtxt("rdf.dat", comments="#")

    assert _float64_sha256(result) == RDF_FLOAT64_SHA256
