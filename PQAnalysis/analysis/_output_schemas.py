"""
Scientific column schemas and xmgrace plot presets for analysis output.
"""

from dataclasses import dataclass



@dataclass(frozen=True)
class AnalysisColumn:

    """Metadata for one analysis-table column."""

    field: str
    symbol: str
    unit: str



@dataclass(frozen=True)
class AnalysisPlotPreset:

    """
    Default xmgrace projection for an analysis table.

    ``x_field=None`` uses a one-based row index as the x axis.
    """

    x_field: str | None
    y_fields: tuple[str, ...]
    x_label: str
    y_label: str
    legends: tuple[str, ...]



@dataclass(frozen=True)
class AnalysisSchema:

    """Title, columns and optional plot defaults for an analysis table."""

    title: str
    columns: tuple[AnalysisColumn, ...]
    plot: AnalysisPlotPreset | None = None

    @property
    def fields(self) -> tuple[str, ...]:
        """Return the stable field identifiers in column order."""
        return tuple(column.field for column in self.columns)

    @property
    def header_columns(self) -> tuple[tuple[str, str, str], ...]:
        """Return columns in the tuple form used by the native header."""
        return tuple(
            (column.field, column.symbol, column.unit)
            for column in self.columns
        )



RDF_SCHEMA = AnalysisSchema(
    title="Radial distribution function",
    columns=(
        AnalysisColumn("r_i", "rᵢ", "Å"),
        AnalysisColumn("g_r_i", "g(rᵢ)", "1"),
        AnalysisColumn("N_r_i", "N(rᵢ)", "1"),
        AnalysisColumn("g_r_i_dV_i", "g(rᵢ)ΔVᵢ", "Å³"),
        AnalysisColumn("H_i_minus_E_i", "Hᵢ−Eᵢ", "pairs"),
    ),
    plot=AnalysisPlotPreset(
        x_field="r_i",
        y_fields=("g_r_i", ),
        x_label=r"r [\cE\C]",
        y_label="g(r)",
        legends=("g(r)", ),
    ),
)

MSD_SCHEMA = AnalysisSchema(
    title="Mean squared displacement",
    columns=(
        AnalysisColumn("lag", "k", "frames"),
        AnalysisColumn("msd_x", "⟨Δx(k)²⟩", "Å²"),
        AnalysisColumn("msd_y", "⟨Δy(k)²⟩", "Å²"),
        AnalysisColumn("msd_z", "⟨Δz(k)²⟩", "Å²"),
    ),
    plot=AnalysisPlotPreset(
        x_field="lag",
        y_fields=("msd_x", "msd_y", "msd_z"),
        x_label="Lag [frames]",
        y_label=r"MSD [\cE\C\S2\N]",
        legends=("x", "y", "z"),
    ),
)

VACF_SCHEMA = AnalysisSchema(
    title="Normalized correlation function",
    columns=(
        AnalysisColumn("time", "t", "ps"),
        AnalysisColumn("normalized_correlation", "C(t)∕C(0)", "1"),
    ),
    plot=AnalysisPlotPreset(
        x_field="time",
        y_fields=("normalized_correlation", ),
        x_label="t [ps]",
        y_label="C(t) / C(0)",
        legends=("C(t) / C(0)", ),
    ),
)

VACF_SPECTRUM_SCHEMA = AnalysisSchema(
    title="VACF spectrum",
    columns=(
        AnalysisColumn("wavenumber", "ν̃", "cm⁻¹"),
        AnalysisColumn("amplitude", "|Ĉ(ν̃)|", "arbitrary"),
    ),
    plot=AnalysisPlotPreset(
        x_field="wavenumber",
        y_fields=("amplitude", ),
        x_label=r"Wavenumber [cm\S-1\N]",
        y_label="Amplitude [a.u.]",
        legends=("Amplitude", ),
    ),
)

VACF_WINDOWED_SCHEMA = AnalysisSchema(
    title="Windowed correlation function",
    columns=(
        AnalysisColumn("time", "t", "ps"),
        AnalysisColumn("windowed_correlation", "C(t)w(t)", "1"),
    ),
    plot=AnalysisPlotPreset(
        x_field="time",
        y_fields=("windowed_correlation", ),
        x_label="t [ps]",
        y_label="C(t)w(t)",
        legends=("Windowed correlation", ),
    ),
)

BROADENED_SPECTRUM_SCHEMA = AnalysisSchema(
    title="Broadened spectrum",
    columns=(
        AnalysisColumn("wavenumber", "ν̃", "cm⁻¹"),
        AnalysisColumn("intensity", "I(ν̃)", "input-dependent"),
    ),
    plot=AnalysisPlotPreset(
        x_field="wavenumber",
        y_fields=("intensity", ),
        x_label=r"Wavenumber [cm\S-1\N]",
        y_label="Intensity",
        legends=("Intensity", ),
    ),
)

MOMENTUM_SCHEMA = AnalysisSchema(
    title="Total linear momentum",
    columns=(
        AnalysisColumn("frame", "n", "1"),
        AnalysisColumn("scaled_momentum_norm", "s‖P(n)‖", "scale-dependent"),
    ),
    plot=AnalysisPlotPreset(
        x_field="frame",
        y_fields=("scaled_momentum_norm", ),
        x_label="Frame",
        y_label="Scaled momentum norm",
        legends=("Total momentum", ),
    ),
)

VIBRATIONAL_SCHEMA = AnalysisSchema(
    title="Vibrational analysis",
    columns=(
        AnalysisColumn("wavenumber", "ν̃ⱼ", "cm⁻¹"),
        AnalysisColumn("force_constant", "kⱼ", "mdyn·Å⁻¹"),
        AnalysisColumn("reduced_mass", "μⱼ", "amu"),
    ),
    plot=AnalysisPlotPreset(
        x_field=None,
        y_fields=("wavenumber", ),
        x_label="Mode",
        y_label=r"Wavenumber [cm\S-1\N]",
        legends=("Wavenumber", ),
    ),
)

VIBRATIONAL_IR_SCHEMA = AnalysisSchema(
    title="Vibrational analysis",
    columns=(
        AnalysisColumn("wavenumber", "ν̃ⱼ", "cm⁻¹"),
        AnalysisColumn("ir_intensity", "Iⱼᴵᴿ", "km·mol⁻¹"),
        AnalysisColumn("force_constant", "kⱼ", "mdyn·Å⁻¹"),
        AnalysisColumn("reduced_mass", "μⱼ", "amu"),
    ),
    plot=AnalysisPlotPreset(
        x_field="wavenumber",
        y_fields=("ir_intensity", ),
        x_label=r"Wavenumber [cm\S-1\N]",
        y_label=r"IR intensity [km mol\S-1\N]",
        legends=("IR intensity", ),
    ),
)

KNOWN_SCHEMAS = {
    schema.fields: schema
    for schema in (
        RDF_SCHEMA,
        MSD_SCHEMA,
        VACF_SCHEMA,
        VACF_SPECTRUM_SCHEMA,
        VACF_WINDOWED_SCHEMA,
        BROADENED_SPECTRUM_SCHEMA,
        MOMENTUM_SCHEMA,
        VIBRATIONAL_SCHEMA,
        VIBRATIONAL_IR_SCHEMA, )
}



def normal_mode_schema(n_modes: int) -> AnalysisSchema:
    """
    Build the dynamic schema for a normal-mode matrix.

    Parameters
    ----------
    n_modes : int
        Number of normal-mode columns.

    Returns
    -------
    AnalysisSchema
        Matrix schema with one field per mode.
    """
    columns = tuple(
        AnalysisColumn(f"mode_{index}", f"e(α,{index})", "1")
        for index in range(1, n_modes + 1)
    )
    fields = tuple(column.field for column in columns)

    return AnalysisSchema(
        title="Normal-mode matrix",
        columns=columns,
        plot=AnalysisPlotPreset(
            x_field=None,
            y_fields=fields,
            x_label="Cartesian component index",
            y_label="Normalized mode component",
            legends=tuple(f"Mode {index}" for index in range(1, n_modes + 1)),
        ),
    )
