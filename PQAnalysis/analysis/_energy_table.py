"""
Column metadata for PQ and QMCFC energy files as analysis tables.
"""

import re
from dataclasses import dataclass
from numbers import Real
from pathlib import Path

from PQAnalysis.physical_data import Energy
from PQAnalysis.traj import MDEngineFormat

from ._output_schemas import AnalysisColumn

ENERGY_FILE_SUFFIX = ".en"
_TIME_FIELDS = ("SIMULATION-TIME", "SIMULATION_TIME", "STEP")
_DIMENSIONLESS_UNITS = (None, "", "-")



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



def energy_plot_fields(
    fields: tuple[str, ...],
) -> tuple[str | None, tuple[str, ...]]:
    """
    Return the default x field (simulation time) and y fields.

    The y fields are all other columns except the step counter.
    """
    x_field = next((name for name in _TIME_FIELDS if name in fields), None)
    if x_field is None and len(fields) > 1:
        x_field = fields[0]

    y_fields = tuple(
        field for field in fields
        if field != x_field and field not in _TIME_FIELDS
    )
    return x_field, y_fields or tuple(
        field for field in fields if field != x_field
    )
