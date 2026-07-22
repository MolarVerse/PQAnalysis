"""
Parity tests of the VACF analysis against the legacy FreqCalc,
Fluxfreqcalc and ft.f tools.

The reference data in tests/data/vacf was generated with the original
(recompiled) legacy binaries from a synthetic 8-atom QMCFC velocity
trajectory (superposed cosines plus noise, split over two files) with
window = 100, gap = 5, time_step = 0.002 ps and target_atoms = 1-8,
and spectra with ftsize = 256 for all window functions with
windowparam = 20.0, winon = 0.02 and winoff = 0.15.

The tolerances account for the float32 parsing of PQAnalysis (the
legacy tools parse in double precision) and for the printing precision
of the reference files.
"""

import numpy as np
import pytest

from PQAnalysis.analysis.vacf import (
    VACF,
    read_static_charges,
    vacf,
    vacf_spectrum,
)
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.traj import TrajectoryFormat

from .. import pytestmark  # pylint: disable=unused-import

VEL_FILES = ["traj_1.vel", "traj_2.vel"]
CHARGE_FILES = ["traj_1.chrg", "traj_2.chrg"]

WINDOW_SIZE = 100
GAP = 5
TIME_STEP = 0.002
FTSIZE = 256

SPECTRUM_KWARGS = {
    "none": {},
    "exponential": {
        "window_function": "exponential",
        "window_param": 20.0,
        "window_start": 0.02,
        "window_stop": 0.15,
    },
    "hann": {
        "window_function": "hann",
        "window_start": 0.02,
        "window_stop": 0.15,
    },
    "blackman": {
        "window_function": "blackman",
        "window_start": 0.02,
        "window_stop": 0.15,
    },
}



def _run_vacf(**kwargs):
    """
    Runs a VACF analysis on the reference velocity trajectory.
    """
    reader = TrajectoryReader(VEL_FILES, md_format="qmcfc")

    analysis = VACF(
        traj=reader,
        window_size=WINDOW_SIZE,
        time_step=TIME_STEP,
        gap=GAP,
        **kwargs,
    )

    return analysis.run()



class TestVACFLegacyParity:

    """
    Parity tests against the legacy reference data.
    """

    @pytest.mark.parametrize("example_dir", ["vacf"], indirect=False)
    def test_freqcalc_parity(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        The VACF matches the legacy FreqCalc reference.
        """
        time, correlation = _run_vacf()

        reference = np.loadtxt("vacf_ref.dat")

        assert np.allclose(time, reference[:, 0], atol=1e-6)
        assert np.allclose(
            correlation,
            reference[:, 1],
            rtol=1e-5,
            atol=1e-8,
        )

    @pytest.mark.parametrize("example_dir", ["vacf"], indirect=False)
    def test_fluxfreqcalc_static_parity(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        The charge-flux VACF with static charges matches the legacy
        Fluxfreqcalc reference.
        """
        charges = read_static_charges("charges.dat", md_format="qmcfc")
        _, correlation = _run_vacf(charges=charges)

        reference = np.loadtxt("flux_static_ref.dat")

        assert np.allclose(
            correlation,
            reference[:, 1],
            rtol=1e-5,
            atol=1e-8,
        )

    @pytest.mark.parametrize("example_dir", ["vacf"], indirect=False)
    def test_fluxfreqcalc_charge_trajectory_parity(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        The charge-flux VACF with a charge trajectory matches the
        legacy Fluxfreqcalc reference.
        """
        charge_reader = TrajectoryReader(
            CHARGE_FILES,
            traj_format=TrajectoryFormat.CHARGE,
            md_format="qmcfc",
        )
        _, correlation = _run_vacf(charge_traj=charge_reader)

        reference = np.loadtxt("flux_traj_ref.dat")

        assert np.allclose(
            correlation,
            reference[:, 1],
            rtol=1e-5,
            atol=1e-8,
        )

    @pytest.mark.parametrize(
        "window_function",
        ["none", "exponential", "hann", "blackman"],
    )
    @pytest.mark.parametrize("example_dir", ["vacf"], indirect=False)
    def test_ft_parity(self, test_with_data_dir, window_function):  # pylint: disable=unused-argument
        """
        The spectrum of the reference VACF matches the legacy ft.f
        reference for all window functions.
        """
        reference_vacf = np.loadtxt("vacf_ref.dat")

        wavenumbers, amplitudes, windowed = vacf_spectrum(
            reference_vacf[:, 0],
            reference_vacf[:, 1],
            ftsize=FTSIZE,
            **SPECTRUM_KWARGS[window_function],
        )

        reference = np.loadtxt(f"spectrum_{window_function}_ref.dat")

        assert np.allclose(wavenumbers, reference[:, 0], atol=1e-6)
        assert np.allclose(
            amplitudes,
            reference[:, 1],
            rtol=1e-5,
            atol=1e-9,
        )

        if window_function != "none":
            windowed_reference = np.loadtxt(
                f"windowed_{window_function}_ref.dat"
            )

            assert np.allclose(
                windowed,
                windowed_reference[:, 1],
                rtol=1e-5,
                atol=1e-9,
            )

    @pytest.mark.parametrize("example_dir", ["vacf"], indirect=False)
    def test_api_end_to_end_parity(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        The input file based api reproduces the legacy VACF, spectrum
        and windowed output files.
        """
        with open("vacf.in", "w", encoding="utf-8") as file:
            file.write(
                "traj_files = [traj_1.vel, traj_2.vel]\n"
                "target_selection = all\n"
                "out_file = vacf_out.dat\n"
                "time_step = 0.002\n"
                "window = 100\n"
                "gap = 5\n"
                "spectrum_file = spectrum_out.dat\n"
                "ftsize = 256\n"
                "window_function = hann\n"
                "window_start = 0.02\n"
                "window_stop = 0.15\n"
                "windowed_out_file = windowed_out.dat\n"
                "log_file = vacf.log\n"
            )

        vacf("vacf.in", md_format="qmcfc")

        result = np.loadtxt("vacf_out.dat")
        reference = np.loadtxt("vacf_ref.dat")

        assert np.allclose(result, reference, rtol=1e-5, atol=1e-7)

        spectrum_result = np.loadtxt("spectrum_out.dat")
        spectrum_reference = np.loadtxt("spectrum_hann_ref.dat")

        # the spectrum of the api is calculated from the full
        # precision VACF, while the legacy ft.f tool reads the VACF
        # rounded to 8 decimals - hence the absolute tolerance
        assert np.allclose(
            spectrum_result,
            spectrum_reference,
            rtol=1e-4,
            atol=1e-6,
        )

        windowed_result = np.loadtxt("windowed_out.dat")
        windowed_reference = np.loadtxt("windowed_hann_ref.dat")

        assert np.allclose(
            windowed_result,
            windowed_reference,
            rtol=1e-4,
            atol=1e-7,
        )

    @pytest.mark.parametrize("example_dir", ["vacf"], indirect=False)
    def test_api_flux_end_to_end_parity(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        The input file based api reproduces the legacy charge-flux
        references for both charge sources.
        """
        with open("flux_static.in", "w", encoding="utf-8") as file:
            file.write(
                "traj_files = [traj_1.vel, traj_2.vel]\n"
                "target_selection = all\n"
                "out_file = flux_static_out.dat\n"
                "time_step = 0.002\n"
                "window = 100\n"
                "gap = 5\n"
                "charge_file = charges.dat\n"
            )

        with open("flux_traj.in", "w", encoding="utf-8") as file:
            file.write(
                "traj_files = [traj_1.vel, traj_2.vel]\n"
                "target_selection = all\n"
                "out_file = flux_traj_out.dat\n"
                "time_step = 0.002\n"
                "window = 100\n"
                "gap = 5\n"
                "charge_files = [traj_1.chrg, traj_2.chrg]\n"
            )

        vacf("flux_static.in", md_format="qmcfc")
        vacf("flux_traj.in", md_format="qmcfc")

        static_result = np.loadtxt("flux_static_out.dat")
        static_reference = np.loadtxt("flux_static_ref.dat")

        assert np.allclose(
            static_result,
            static_reference,
            rtol=1e-5,
            atol=1e-7,
        )

        traj_result = np.loadtxt("flux_traj_out.dat")
        traj_reference = np.loadtxt("flux_traj_ref.dat")

        assert np.allclose(
            traj_result,
            traj_reference,
            rtol=1e-5,
            atol=1e-7,
        )
