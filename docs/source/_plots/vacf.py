"""Analytical damped VACF and its PQAnalysis spectrum."""

import matplotlib.pyplot as plt
import numpy as np

from PQAnalysis.analysis.vacf.spectrum import vacf_spectrum

from _style import COLORS, apply_style


apply_style((6.2, 5.5))

time = np.arange(0.0, 0.3005, 0.0005)
correlation = (
    0.85 * np.exp(-(time / 0.075) ** 2)
    * np.cos(2.0 * np.pi * 9.0 * time)
    + 0.15 * np.exp(-(time / 0.035) ** 2)
    * np.cos(2.0 * np.pi * 42.0 * time)
)

wavenumbers, amplitudes, _ = vacf_spectrum(
    time,
    correlation,
    ftsize=4096,
    window_function="hann",
    window_stop=float(time[-1]),
)
display_range = wavenumbers <= 4000.0
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
)
correlation_axis.axhline(
    0.0,
    color=COLORS["muted"],
    linestyle=":",
    linewidth=1.0,
)
correlation_axis.set_title(
    "(a) Normalized damped VACF",
    loc="left",
    fontsize=9.5,
    fontweight="bold",
    pad=8,
)
correlation_axis.set_xlabel("Lag time, t / ps")
correlation_axis.set_ylabel("Cᵥᵥ(t)")
correlation_axis.set_xlim(time[0], time[-1])
correlation_axis.set_ylim(-0.65, 1.05)

spectrum_axis.plot(
    wavenumbers,
    amplitudes,
    color=COLORS["orange"],
)
spectrum_axis.set_title(
    "(b) Hann-window spectrum",
    loc="left",
    fontsize=9.5,
    fontweight="bold",
    pad=8,
)
spectrum_axis.set_xlabel("Wavenumber, ν̃ / cm⁻¹")
spectrum_axis.set_ylabel("Relative amplitude")
spectrum_axis.set_xlim(0.0, 4000.0)
spectrum_axis.set_ylim(0.0, 1.05)

figure.tight_layout(h_pad=1.4)
plt.show()
