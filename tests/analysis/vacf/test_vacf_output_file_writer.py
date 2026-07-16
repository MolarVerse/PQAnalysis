"""
Tests for the VACF output file writers.
"""

import numpy as np

from PQAnalysis.analysis.vacf import (
    VACF,
    VACFDataWriter,
    VACFLogWriter,
    VACFSpectrumDataWriter,
    VACFWindowedDataWriter,
)
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.traj import Trajectory

from .. import pytestmark  # pylint: disable=unused-import

# pylint: disable=protected-access



def _read_lines(filename):
    """
    Reads all lines of a file without trailing newlines.
    """
    with open(filename, "r", encoding="utf-8") as file:
        return [line.rstrip("\n") for line in file]



def _make_vacf():
    """
    Builds a small VACF analysis object.
    """
    rng = np.random.default_rng(42)
    atoms = [Atom("O")]

    systems = [
        AtomicSystem(atoms=atoms, vel=rng.standard_normal((1, 3)))
        for _ in range(20)
    ]

    return VACF(Trajectory(systems), window_size=5, time_step=0.1)



class TestVACFDataWriter:

    """
    Tests for the VACFDataWriter class.
    """

    def test_write_legacy_format(self, tmpdir):  # pylint: disable=unused-argument
        """
        The data writer replicates the legacy FreqCalc output format.
        """
        writer = VACFDataWriter("vacf.dat")
        writer.write(
            (np.array([0.0, 0.002]), np.array([1.0, 0.81100607]))
        )

        assert _read_lines("vacf.dat") == [
            "  0.000000      1.00000000",
            "  0.002000      0.81100607",
        ]



class TestVACFSpectrumDataWriter:

    """
    Tests for the VACFSpectrumDataWriter class.
    """

    def test_write_legacy_format(self, tmpdir):  # pylint: disable=unused-argument
        """
        The spectrum writer replicates the legacy ft.f output format.
        """
        writer = VACFSpectrumDataWriter("spectrum.dat")
        writer.write(
            (
                np.array([32.7023623, 65.4047245]),
                np.array([0.0904911818, 0.4871028398]),
            )
        )

        assert _read_lines("spectrum.dat") == [
            "   32.7023623    0.0904911818",
            "   65.4047245    0.4871028398",
        ]



class TestVACFWindowedDataWriter:

    """
    Tests for the VACFWindowedDataWriter class.
    """

    def test_write_legacy_format(self, tmpdir):  # pylint: disable=unused-argument
        """
        The windowed writer replicates the legacy ft.f windowfile
        format.
        """
        writer = VACFWindowedDataWriter("windowed.dat")
        writer.write(
            (np.array([0.0, 0.002]), np.array([1.0, 0.81100607]))
        )

        assert _read_lines("windowed.dat") == [
            "   0.0000    1.0000000000",
            "   0.0020    0.8110060700",
        ]



class TestVACFLogWriter:

    """
    Tests for the VACFLogWriter class.
    """

    def test_write_before_and_after_run(self, tmpdir):  # pylint: disable=unused-argument
        """
        The log writer writes the setup parameters before the run and
        the number of origins and elapsed time after the run.
        """
        vacf = _make_vacf()

        writer = VACFLogWriter("vacf.log")
        writer.write_before_run(vacf)

        vacf.run()

        writer.write_after_run(vacf)

        log_content = "\n".join(_read_lines("vacf.log"))

        assert "VACF calculation:" in log_content
        assert "Window size (frames): 5" in log_content
        assert "Origin gap (frames):  1" in log_content
        assert "Time step:            0.1 ps" in log_content
        assert "Method:               direct" in log_content
        assert "Number of frames: 20" in log_content
        assert "Number of atoms:  1" in log_content
        assert f"Number of origins: {vacf.n_origins}" in log_content
        assert f"Elapsed time: {vacf.elapsed_time} s" in log_content

    def test_flux_header(self, tmpdir):  # pylint: disable=unused-argument
        """
        In the charge-flux mode the log header names the charge-flux
        auto-correlation.
        """
        rng = np.random.default_rng(43)
        atoms = [Atom("O")]

        systems = [
            AtomicSystem(atoms=atoms, vel=rng.standard_normal((1, 3)))
            for _ in range(20)
        ]

        vacf = VACF(
            Trajectory(systems),
            window_size=5,
            time_step=0.1,
            charges=np.array([-0.8]),
        )

        writer = VACFLogWriter("flux.log")
        writer.write_before_run(vacf)

        log_content = "\n".join(_read_lines("flux.log"))

        assert "Charge-flux auto-correlation calculation:" in log_content
