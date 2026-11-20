"""
Column metadata for PQ and QMCFC energy files as analysis tables.
"""

import re
from dataclasses import dataclass
from numbers import Real
from pathlib import Path

import numpy as np

from PQAnalysis.physical_data import Energy
from PQAnalysis.traj import MDEngineFormat

from ._output_schemas import (
    AnalysisColumn,
    AnalysisPlotGraph,
    AnalysisPlotPreset,
    AnalysisSchema,
)

ENERGY_FILE_SUFFIX = ".en"
_TIME_FIELDS = ("SIMULATION-TIME", "SIMULATION_TIME", "STEP")
_DIMENSIONLESS_UNITS = (None, "", "-")
MAX_DEFAULT_GRAPHS = 4
_ENERGY_UNITS = ("kcal/mol", "kj/mol", "ev", "hartree", "ha", "kcal", "kj")
_ENERGY_PREFIXES = ("E(", "E_")
_QUANTITY_BY_UNIT = {
    "K": "Temperature",
    "bar": "Pressure",
    "A^3": "Volume",
    "g/cm^3": "Density",
    "s": "Time",
    "ps": "Time",
}



@dataclass(frozen=True)
class EnergyFileOptions:

    """
    How to read an energy-file input of ``convert_analysis_output``.

    Attributes
    ----------
    as_energy_file : bool
        Read the input as an energy file whatever its suffix. Inputs
        ending in ``.en`` are always read as energy files.
    info_file : str | None
        Info file with the column names and units, by default the
        energy filename with the ``.info`` suffix if it exists.
    engine_format : MDEngineFormat | str
        Engine that wrote the file, by default PQ.
    timestep_fs : Real | None
        PQ integration timestep in fs, used to report the simulation
        time in ps instead of steps.
    """

    as_energy_file: bool = False
    info_file: str | None = None
    engine_format: MDEngineFormat | str = "PQ"
    timestep_fs: Real | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "engine_format", MDEngineFormat(self.engine_format)
        )

    def reads_energy(self, filename: str) -> bool:
        """Return whether ``filename`` is read as an energy file."""
        return self.as_energy_file or is_energy_file(filename)

    @property
    def has_file_options(self) -> bool:
        """Return whether an info file or a timestep was given."""
        return self.info_file is not None or self.timestep_fs is not None



def is_energy_file(filename: str) -> bool:
    """
    Return whether a filename has the PQ/QMCFC energy-file suffix.
    """
    return Path(filename).suffix.lower() == ENERGY_FILE_SUFFIX



def energy_columns(energy: Energy) -> tuple[AnalysisColumn, ...]:
    """
    Return one analysis column per energy-file column, in file order.

    Names come from the energy info; whitespace becomes ``_`` and
    dimensionless units (``-`` or none) become ``1``. Without info the
    fields are ``column_1`` ... ``column_n``.
    """
    n_columns = len(energy.data)
    if not energy.info_given:
        return tuple(
            AnalysisColumn(f"column_{index}", f"column_{index}", "1")
            for index in range(1, n_columns + 1)
        )

    units = energy.units if energy.units_given else {}
    columns = []
    for name in sorted(energy.info, key=energy.info.get):
        field = re.sub(r"\s+", "_", name.strip())
        unit = units.get(name)
        unit = "1" if unit in _DIMENSIONLESS_UNITS else str(unit)
        columns.append(AnalysisColumn(field, field, unit))
    return tuple(columns)



def default_energy_plot(
    schema: AnalysisSchema,
    data: np.ndarray,
) -> AnalysisPlotPreset:
    """
    Return the default plot: one stacked graph per energy column.

    The graphs show the first ``MAX_DEFAULT_GRAPHS`` energy columns
    (energy unit, or an ``E(`` name if the file has no units) that are
    not zero everywhere, each on its own y axis because total and
    kinetic energies differ by orders of magnitude. Without such
    columns, all columns except the time and step counter, grouped by
    unit.
    """
    x_field, y_fields, by_energy = _default_plot_fields(schema.columns, data)
    graphs = tuple((field, ) for field in y_fields) if by_energy else None
    return energy_plot_preset(schema, x_field, y_fields, graphs)



def _default_plot_fields(
    columns: tuple[AnalysisColumn, ...],
    data: np.ndarray,
) -> tuple[str | None, tuple[str, ...], bool]:
    fields = tuple(column.field for column in columns)
    x_field = next((name for name in _TIME_FIELDS if name in fields), None)
    if x_field is None and len(fields) > 1:
        x_field = fields[0]

    energies = tuple(
        column.field
        for index, column in enumerate(columns)
        if column.field != x_field and _is_energy(column) and
        np.any(data[:, index] != 0)
    )
    if energies:
        return x_field, energies[:MAX_DEFAULT_GRAPHS], True

    y_fields = tuple(
        field for field in fields
        if field != x_field and field not in _TIME_FIELDS
    )
    return x_field, y_fields or tuple(
        field for field in fields if field != x_field
    ), False



def energy_plot_preset(
    schema: AnalysisSchema,
    x_field: str | None,
    y_fields: tuple[str, ...],
    graphs: tuple[tuple[str, ...], ...] | None = None,
) -> AnalysisPlotPreset:
    """
    Return a plot preset with stacked graphs for the y fields.

    Without explicit ``graphs`` there is one graph per unit, in order of
    first appearance in ``y_fields``; without units, energies and other
    quantities get separate graphs.
    """
    by_field = {column.field: column for column in schema.columns}
    if graphs is None:
        groups: dict[tuple[str, bool], list[str]] = {}
        for field in y_fields:
            column = by_field[field]
            key = (column.unit, _is_energy(column))
            groups.setdefault(key, []).append(field)
        graphs = tuple(tuple(fields) for fields in groups.values())

    plot_graphs = tuple(
        AnalysisPlotGraph(
            y_label=_graph_label([by_field[field] for field in fields]),
            y_fields=fields,
        ) for fields in graphs
    )
    if x_field is None:
        x_label = "Row"
    elif x_field in _TIME_FIELDS:
        x_label = _unit_suffix("Simulation time", by_field[x_field].unit)
    else:
        x_label = energy_legend(by_field[x_field])

    return AnalysisPlotPreset(
        x_field=x_field,
        y_fields=y_fields,
        x_label=x_label,
        y_label=plot_graphs[0].y_label,
        legends=y_fields,
        graphs=plot_graphs,
    )



def energy_legend(column: AnalysisColumn) -> str:
    """Return ``NAME [unit]`` in Grace notation, the unit if not ``1``."""
    return _unit_suffix(column.field, column.unit)



def _is_energy(column: AnalysisColumn) -> bool:
    if column.unit in _DIMENSIONLESS_UNITS + ("1", ):
        return column.field.startswith(_ENERGY_PREFIXES)
    return column.unit.lower() in _ENERGY_UNITS



def _graph_label(columns: list[AnalysisColumn]) -> str:
    if len(columns) == 1:
        return energy_legend(columns[0])
    if all(_is_energy(column) for column in columns):
        return _unit_suffix("Energy", columns[0].unit)
    quantity = _QUANTITY_BY_UNIT.get(columns[0].unit, "Value")
    return _unit_suffix(quantity, columns[0].unit)



def _unit_suffix(label: str, unit: str) -> str:
    if unit in _DIMENSIONLESS_UNITS + ("1", ):
        return label
    return f"{label} [{grace_unit(unit)}]"



def grace_unit(unit: str) -> str:
    """Return a unit in Grace notation: ``A`` as Angstrom, ``^n`` raised."""
    unit = re.sub(r"(?<![A-Za-z])A(?=\^|$)", r"\\cE\\C", unit)
    return re.sub(r"\^(-?\d+)", r"\\S\1\\N", unit)
