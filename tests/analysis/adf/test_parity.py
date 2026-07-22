"""
Parity tests of the PQAnalysis ADF analysis against the native C
reference ``adf_ref`` (a double-precision port of the legacy thh_tools
ADF algorithm, source in ``tests/data/adf/adf_ref.c``).

The reference histograms committed in ``tests/data/adf`` were produced
by ``adf_ref`` on ``traj.xyz`` (an orthorhombic 8-frame trajectory).
PQAnalysis parses positions in float32 and images them with the general
triclinic minimum-image convention, whereas the reference works in
double precision with the legacy orthorhombic ``round`` imaging, so a
handful of triplets sitting exactly on an angle bin edge are allowed to
differ. On this trajectory the integer histograms match exactly.
"""

import numpy as np
import pytest

from PQAnalysis.analysis import ADF
from PQAnalysis.io import TrajectoryReader

from .. import pytestmark  # pylint: disable=unused-import


def _reference_counts(filename):
    reference = np.loadtxt(filename)
    return reference[:, 0], reference[:, 1].astype(np.int64)



class TestADFParity:

    @pytest.mark.parametrize("example_dir", ["adf"], indirect=False)
    def test_parity_no_gate(self, test_with_data_dir):
        centers, ref_counts = _reference_counts("adf_ref.out")

        adf = ADF(TrajectoryReader("traj.xyz"), "O", "H", n_angle_bins=180)
        mids, _norm, counts, _sinc = adf.run()
        counts = counts.astype(np.int64)

        assert np.allclose(mids, centers)
        # same total number of counted (ordered) triplets
        assert counts.sum() == ref_counts.sum() == 30400
        # the integer histograms match exactly on this trajectory
        assert np.array_equal(counts, ref_counts)

    @pytest.mark.parametrize("example_dir", ["adf"], indirect=False)
    def test_parity_gated(self, test_with_data_dir):
        centers, ref_counts = _reference_counts("adf_ref_gated.out")

        adf = ADF(
            TrajectoryReader("traj.xyz"),
            "O",
            "H",
            delta_angle=2.0,
            r_min_1=0.8,
            r_max_1=3.0,
            r_min_2=0.8,
            r_max_2=3.0,
        )
        mids, _norm, counts, _sinc = adf.run()
        counts = counts.astype(np.int64)

        assert np.allclose(mids, centers)
        assert counts.sum() == ref_counts.sum()
        assert np.array_equal(counts, ref_counts)

    @pytest.mark.parametrize("example_dir", ["adf"], indirect=False)
    def test_parity_bins_match_within_tolerance(self, test_with_data_dir):
        # quantified parity: at most a negligible fraction of triplets
        # may cross a bin edge due to the float32-vs-float64 difference
        _centers, ref_counts = _reference_counts("adf_ref.out")

        adf = ADF(TrajectoryReader("traj.xyz"), "O", "H", n_angle_bins=180)
        _mids, _norm, counts, _sinc = adf.run()
        counts = counts.astype(np.int64)

        total = ref_counts.sum()
        sum_abs_diff = int(np.abs(counts - ref_counts).sum())

        # fewer than 0.01 % of the counted triplets may land in a
        # neighbouring bin (here: zero)
        assert sum_abs_diff <= total * 1e-4
