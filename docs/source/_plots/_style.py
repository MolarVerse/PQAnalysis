"""Shared Matplotlib style for the scientific documentation figures."""

from pathlib import Path

import matplotlib as mpl


PROJECT_ROOT = Path(__file__).resolve().parents[3]

COLORS = {
    "blue": "#176c8c",
    "green": "#008f72",
    "orange": "#c7521c",
    "magenta": "#a84d84",
    "ink": "#202428",
    "muted": "#66717a",
    "grid": "#d7dde1",
    "shell": "#dcecf2",
}


def apply_style(figsize: tuple[float, float]) -> None:
    """Apply a restrained, colorblind-safe style to one figure."""

    mpl.rcParams.update({
        "figure.figsize": figsize,
        "figure.dpi": 120,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
        "font.size": 9.5,
        "axes.labelsize": 10,
        "axes.labelcolor": COLORS["ink"],
        "axes.edgecolor": COLORS["muted"],
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.axisbelow": True,
        "axes.grid": True,
        "grid.color": COLORS["grid"],
        "grid.linewidth": 0.7,
        "grid.alpha": 0.8,
        "xtick.color": COLORS["ink"],
        "ytick.color": COLORS["ink"],
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "lines.linewidth": 1.8,
    })
