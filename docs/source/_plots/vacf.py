"""VACF and Hann-window spectrum from the validation fixture."""

import matplotlib.pyplot as plt
import numpy as np

from _style import COLORS, PROJECT_ROOT, apply_style


apply_style((6.2, 5.5))

correlation = np.loadtxt(PROJECT_ROOT / "tests/data/vacf/vacf_ref.dat")
spectrum = np.loadtxt(
    PROJECT_ROOT / "tests/data/vacf/spectrum_hann_ref.dat"
)
spectrum = spectrum[spectrum[:, 0] <= 4000.0]

figure, (correlation_axis, spectrum_axis) = plt.subplots(
    2,
    1,
    gridspec_kw={"height_ratios": (1.25, 1.0)},
)

correlation_axis.plot(
    correlation[:, 0],
    correlation[:, 1],
    color=COLORS["blue"],
)
correlation_axis.axhline(
    0.0,
    color=COLORS["muted"],
    linestyle=":",
    linewidth=1.0,
)
correlation_axis.set_title(
    "(a) Normalized velocity autocorrelation",
    loc="left",
    fontsize=9.5,
    fontweight="bold",
    pad=8,
)
correlation_axis.set_xlabel("Lag time, t / ps")
correlation_axis.set_ylabel("Cᵥᵥ(t)")
correlation_axis.set_xlim(correlation[0, 0], correlation[-1, 0])
correlation_axis.set_ylim(-1.05, 1.05)

spectrum_axis.plot(
    spectrum[:, 0],
    spectrum[:, 1],
    color=COLORS["orange"],
)
spectrum_axis.set_title(
    "(b) Hann-window cosine-transform spectrum",
    loc="left",
    fontsize=9.5,
    fontweight="bold",
    pad=8,
)
spectrum_axis.set_xlabel("Wavenumber, ν̃ / cm⁻¹")
spectrum_axis.set_ylabel("|Ĉ(ν̃)| / a.u.")
spectrum_axis.set_xlim(0.0, 4000.0)
spectrum_axis.set_ylim(bottom=0.0)

figure.tight_layout(h_pad=1.4)
plt.show()
