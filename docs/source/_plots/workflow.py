"""RDF and MSD of the bundled water tutorial fixture."""

import matplotlib.pyplot as plt

from PQAnalysis import config
from PQAnalysis.analysis import MSD, RDF
from PQAnalysis.io import TrajectoryReader

from _style import COLORS, PROJECT_ROOT, apply_style

config.with_progress_bar = False

water = PROJECT_ROOT / "examples" / "water"
xyz = str(water / "trajectory.xyz")

r, g, coordination, _shell, _residual = RDF(
    TrajectoryReader(xyz),
    reference_species="O",
    target_species="H",
    delta_r=0.5,
    r_max=4.0,
).run()

lags, _msd_x, _msd_y, _msd_z, msd_tot = MSD(
    TrajectoryReader(xyz),
    target_species="O",
    window=8,
    gap=2,
    time_step=0.001,
    fit_window=4,
).run()
time = lags * 0.001

apply_style((6.4, 5.4))
figure, (rdf_axis, msd_axis) = plt.subplots(
    2,
    1,
    gridspec_kw={"height_ratios": (1.15, 1.0)},
)

rdf_axis.plot(r, g, color=COLORS["blue"], label=r"$g(r)$")
rdf_axis.set_ylabel(r"$g(r)$")
rdf_axis.set_xlabel(r"$r$ / Å")
rdf_axis.legend(loc="upper right")

msd_axis.plot(time, msd_tot, color=COLORS["green"], label="total MSD")
msd_axis.set_ylabel(r"MSD / Å$^2$")
msd_axis.set_xlabel(r"$t$ / ps")
msd_axis.legend(loc="upper left")

figure.tight_layout()
