"""MSD components and diffusion-fit interval from the validation fixture."""

import matplotlib.pyplot as plt
import numpy as np

from _style import COLORS, PROJECT_ROOT, apply_style


apply_style((7.2, 4.3))

data = np.loadtxt(PROJECT_ROOT / "tests/data/msd/msd_ref_O.dat")
time = data[:, 0] * 0.5
components = data[:, 1:4]
total = np.sum(components, axis=1)

fit_start = len(time) - 20
fit_coefficients = np.polyfit(time[fit_start:], total[fit_start:], 1)
fit = np.polyval(fit_coefficients, time[fit_start:])

figure, axis = plt.subplots()
for values, label, color in zip(
    components.T,
    (r"$\mathrm{MSD}_x$", r"$\mathrm{MSD}_y$", r"$\mathrm{MSD}_z$"),
    (COLORS["blue"], COLORS["green"], COLORS["magenta"]),
):
    axis.plot(time, values, color=color, linewidth=1.35, label=label)

axis.plot(time, total, color=COLORS["ink"], linewidth=2.2, label="total")
axis.axvspan(
    time[fit_start],
    time[-1],
    color=COLORS["shell"],
    label="fit interval",
)
axis.plot(
    time[fit_start:],
    fit,
    color=COLORS["orange"],
    linestyle="--",
    linewidth=1.7,
    label="linear fit",
)
axis.set_xlabel(r"Lag time $t$ / ps")
axis.set_ylabel(r"Mean square displacement / $\mathrm{\AA}^2$")
axis.set_xlim(time[0], time[-1])
axis.set_ylim(bottom=0.0)
axis.legend(ncol=3, loc="upper left")

figure.tight_layout()
plt.show()
