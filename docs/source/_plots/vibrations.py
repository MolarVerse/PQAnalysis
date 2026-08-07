"""IR stick spectrum calculated from the H2O validation fixture."""

import matplotlib.pyplot as plt
import numpy as np

from PQAnalysis.analysis.vibrational.vibrational_analysis import (
    calculate_from_system,
    read_hessian_file,
)
from PQAnalysis.io import MoldescriptorReader, RestartFileReader

from _style import COLORS, PROJECT_ROOT, apply_style


apply_style((7.2, 3.8))

fixture = PROJECT_ROOT / "tests/data/vibrational"
moldescriptor = fixture / "moldescriptor.dat"
system = RestartFileReader(
    str(fixture / "h2o.rst"),
    moldescriptor_filename=str(moldescriptor),
).read()
hessian = read_hessian_file(str(fixture / "hessian.dat"))
charges = np.asarray(
    MoldescriptorReader(str(moldescriptor)).read()[0].partial_charges,
    dtype=float,
)
result = calculate_from_system(
    system,
    hessian,
    atom_charges=charges,
)

internal_modes = result.wavenumbers > 100.0
wavenumbers = result.wavenumbers[internal_modes]
intensities = result.intensities[internal_modes]

figure, axis = plt.subplots()
axis.vlines(
    wavenumbers,
    0.0,
    intensities,
    color=COLORS["blue"],
    linewidth=2.0,
)
axis.scatter(
    wavenumbers,
    intensities,
    color=COLORS["blue"],
    marker="_",
    s=85,
)
for index, (wavenumber, intensity) in enumerate(
    zip(wavenumbers, intensities)
):
    axis.annotate(
        f"{wavenumber:.0f}",
        (wavenumber, intensity),
        xytext=(0, 5 + 10 * (index % 2)),
        textcoords="offset points",
        ha="center",
        va="bottom",
        color=COLORS["ink"],
        fontsize=8,
    )

axis.set_xlabel(r"Wavenumber $\tilde{\nu}$ / $\mathrm{cm}^{-1}$")
axis.set_ylabel(r"IR intensity / $\mathrm{km\ mol}^{-1}$")
axis.set_xlim(0.0, 4200.0)
axis.set_ylim(0.0, max(intensities) * 1.22)

figure.tight_layout()
plt.show()
