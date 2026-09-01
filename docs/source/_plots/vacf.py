"""Analytical two-band VACF and its PQAnalysis spectrum."""

import matplotlib.pyplot as plt
import numpy as np

from PQAnalysis.analysis.vacf.spectrum import vacf_spectrum

from _style import COLORS, apply_style


apply_style((6.2, 5.5))

time = np.arange(0.0, 0.5005, 0.0005)
correlation = (
    0.70 * np.exp(-(time / 0.22) ** 2)
    * np.cos(2.0 * np.pi * 9.0 * time)
    + 0.30 * np.exp(-(time / 0.12) ** 2)
    * np.cos(2.0 * np.pi * 18.0 * time)
)

wavenumbers, amplitudes, windowed_correlation = vacf_spectrum(
    time,
    correlation,
    ftsize=5000,
    window_function="exponential",
    window_param=4.0,
)
display_range = wavenumbers <= 1000.0
wavenumbers = wavenumbers[display_range]
amplitudes = amplitudes[display_range]
amplitudes /= amplitudes.max()

figure, (correlation_axis, spectrum_axis) = plt.subplots(
    2,
    1,
    gridspec_kw={"height_ratios": (1.25, 1.0)},
)

correlation_axis.plot(
    time,
    correlation,
    color=COLORS["blue"],
    label="Unwindowed",
)
correlation_axis.plot(
    time,
    windowed_correlation,
    color=COLORS["green"],
    linestyle="--",
    linewidth=1.6,
    label="Exponential, 4 ps⁻¹",
)
correlation_axis.axhline(
    0.0,
    color=COLORS["muted"],
    linestyle=":",
    linewidth=1.0,
)
correlation_axis.set_title(
    "(a) Normalized VACF",
    loc="left",
    fontsize=9.5,
    fontweight="bold",
    pad=8,
)
correlation_axis.set_xlabel("Lag time, t / ps")
correlation_axis.set_ylabel("Cᵥᵥ(t)")
correlation_axis.set_xlim(time[0], time[-1])
correlation_axis.set_ylim(-0.65, 1.05)
correlation_axis.legend(loc="upper right")

spectrum_axis.plot(
    wavenumbers,
    amplitudes,
    color=COLORS["orange"],
)
spectrum_axis.set_title(
    "(b) Exponential-window spectrum",
    loc="left",
    fontsize=9.5,
    fontweight="bold",
    pad=8,
)
spectrum_axis.set_xlabel("Wavenumber, ν̃ / cm⁻¹")
spectrum_axis.set_ylabel("Relative amplitude")
spectrum_axis.set_xlim(0.0, 1000.0)
spectrum_axis.set_ylim(0.0, 1.05)

figure.tight_layout(h_pad=1.4)
plt.show()
