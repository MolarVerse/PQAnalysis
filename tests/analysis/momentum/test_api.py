"""
Tests for the momentum analysis API and output writer.
"""

import numpy as np
import pytest

from PQAnalysis.analysis.momentum import MomentumDataWriter, check_momentum

from .. import pytestmark  # pylint: disable=unused-import

from .test_momentum import TWO_FRAMES_NORMS



class TestCheckMomentumAPI:

    """
    Tests for the check_momentum API function.
    """

    @pytest.mark.parametrize("example_dir", ["momentum"], indirect=False)
    def test_check_momentum_to_file(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        The API reads a QMCFC velocity trajectory and writes one row
        per frame with the one-based frame index and the scaled
        momentum norm.
        """
        momentum_norms = check_momentum(
            "two_frames.vel",
            output="momentum.dat",
            md_format="qmcfc",
        )

        assert np.allclose(
            momentum_norms,
            np.array(TWO_FRAMES_NORMS) * 1e-15,
            rtol=1e-12,
        )

        with open("momentum.dat", encoding="utf-8") as file:
            lines = file.readlines()

        assert lines == [
            f"1  {momentum_norms[0]:.12e}\n",
            f"2  {momentum_norms[1]:.12e}\n",
        ]

    @pytest.mark.parametrize("example_dir", ["momentum"], indirect=False)
    def test_check_momentum_to_stdout(self, test_with_data_dir, capsys):  # pylint: disable=unused-argument
        """
        Without an output file the momentum norms are printed
        to stdout.
        """
        momentum_norms = check_momentum(
            "two_frames.vel",
            md_format="qmcfc",
            scale=1.0,
        )

        captured = capsys.readouterr()

        assert captured.out == (
            f"1  {momentum_norms[0]:.12e}\n"
            f"2  {momentum_norms[1]:.12e}\n"
        )
        assert np.allclose(momentum_norms, TWO_FRAMES_NORMS, rtol=1e-12)



class TestMomentumDataWriter:

    """
    Tests for the MomentumDataWriter class.
    """

    def test_write(self, tmpdir):  # pylint: disable=unused-argument
        """
        The writer produces one-based frame indices and the momentum
        norms in scientific notation.
        """
        writer = MomentumDataWriter("momentum_out.dat")
        writer.write(np.array([1.871482266947e-14, 0.0]))

        with open("momentum_out.dat", encoding="utf-8") as file:
            lines = file.readlines()

        assert lines == [
            "1  1.871482266947e-14\n",
            "2  0.000000000000e+00\n",
        ]
