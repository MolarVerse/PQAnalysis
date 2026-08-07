"""VACF and Hann-window spectrum from the validation fixture."""

import matplotlib.pyplot as plt
import numpy as np

from _style import COLORS, PROJECT_ROOT, apply_style


apply_style((7.2, 5.1))

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
correlation_axis.set_xlabel(r"Lag time $t$ / ps")
correlation_axis.set_ylabel(r"$C_v(t)$")
correlation_axis.set_xlim(correlation[0, 0], correlation[-1, 0])
correlation_axis.set_ylim(-1.05, 1.05)

spectrum_axis.plot(
    spectrum[:, 0],
    spectrum[:, 1],
    color=COLORS["orange"],
)
spectrum_axis.set_xlabel(r"Wavenumber $\tilde{\nu}$ / $\mathrm{cm}^{-1}$")
spectrum_axis.set_ylabel("Amplitude / a.u.")
spectrum_axis.set_xlim(0.0, 4000.0)
spectrum_axis.set_ylim(bottom=0.0)

figure.tight_layout()
plt.show()
