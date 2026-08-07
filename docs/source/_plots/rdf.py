"""Analytic RDF profile used to explain structural features."""

import matplotlib.pyplot as plt
import numpy as np

from _style import COLORS, apply_style


apply_style((7.2, 5.0))

r = np.linspace(0.02, 8.0, 800)
excluded_volume = 1.0 - np.exp(-(r / 1.65)**8)
structure = (
    1.0
    + 2.2 * np.exp(-0.5 * ((r - 2.80) / 0.22)**2)
    - 0.55 * np.exp(-0.5 * ((r - 3.55) / 0.30)**2)
    + 0.65 * np.exp(-0.5 * ((r - 4.65) / 0.38)**2)
    - 0.18 * np.exp(-0.5 * ((r - 5.55) / 0.45)**2)
)
g_r = np.clip(excluded_volume * structure, 0.0, None)

number_density = 0.0334
coordination_integrand = 4.0 * np.pi * number_density * r**2 * g_r
coordination = np.concatenate((
    [0.0],
    np.cumsum(
        0.5
        * (coordination_integrand[1:] + coordination_integrand[:-1])
        * np.diff(r)
    ),
))

first_minimum = 3.55
figure, (rdf_axis, coordination_axis) = plt.subplots(
    2,
    1,
    sharex=True,
    gridspec_kw={"height_ratios": (2.0, 1.25)},
)

rdf_axis.axvspan(
    0.0,
    first_minimum,
    color=COLORS["shell"],
    label="first coordination shell",
)
rdf_axis.plot(r, g_r, color=COLORS["blue"])
rdf_axis.axhline(1.0, color=COLORS["muted"], linestyle=":", linewidth=1.1)
rdf_axis.axvline(
    first_minimum,
    color=COLORS["orange"],
    linestyle="--",
    linewidth=1.2,
    label="first minimum",
)
rdf_axis.set_ylabel(r"$g(r)$")
rdf_axis.set_ylim(0.0, 3.6)
rdf_axis.legend(loc="upper right")

coordination_axis.plot(r, coordination, color=COLORS["orange"])
coordination_axis.axvline(
    first_minimum,
    color=COLORS["orange"],
    linestyle="--",
    linewidth=1.2,
)
coordination_axis.set_xlabel(r"Distance $r$ / $\mathrm{\AA}$")
coordination_axis.set_ylabel(r"$N(r)$")
coordination_axis.set_xlim(0.0, 8.0)
coordination_axis.set_ylim(bottom=0.0)

figure.tight_layout()
plt.show()
