"""
Read, write and convert self-describing analysis tables.
"""

from __future__ import annotations

import csv
import io
from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace
from enum import Enum
from numbers import Real
from pathlib import Path
from typing import TextIO

import numpy as np

from PQAnalysis.exceptions import PQException
from PQAnalysis.io import BaseWriter
from PQAnalysis.io.formats import FileWritingMode

from ._output_schemas import (
    AnalysisColumn,
    AnalysisPlotPreset,
    AnalysisSchema,
    BROADENED_SPECTRUM_SCHEMA,
    KNOWN_SCHEMAS,
    MOMENTUM_SCHEMA,
    MSD_SCHEMA,
    RDF_SCHEMA,
    VACF_SCHEMA,
    VACF_SPECTRUM_SCHEMA,
    VACF_WINDOWED_SCHEMA,
    VIBRATIONAL_IR_SCHEMA,
    VIBRATIONAL_SCHEMA,
    normal_mode_schema,
)
from ._xvg import (
    XVGDocument,
    XVGFormatError,
    XVGPlotStyle,
    is_xvg,
    read_xvg,
    write_xvg,
)

__all__ = [
    "AnalysisColumn",
    "AnalysisDataWriter",
    "AnalysisOutputError",
    "AnalysisOutputFormat",
    "AnalysisPlotPreset",
    "AnalysisSchema",
    "AnalysisTable",
    "AnalysisTableWriter",
    "BROADENED_SPECTRUM_SCHEMA",
    "MOMENTUM_SCHEMA",
    "MSD_SCHEMA",
    "RDF_SCHEMA",
    "VACF_SCHEMA",
    "VACF_SPECTRUM_SCHEMA",
    "VACF_WINDOWED_SCHEMA",
    "VIBRATIONAL_IR_SCHEMA",
    "VIBRATIONAL_SCHEMA",
    "convert_analysis_output",
    "infer_output_format",
    "normal_mode_schema",
    "read_analysis_table",
    "write_analysis_table",
]

AnalysisTextStream = TextIO | io.StringIO



class AnalysisOutputError(PQException):

    """
    Raised when an analysis table cannot be read, written or converted.
    """



class AnalysisOutputFormat(str, Enum):

    """
    Supported analysis-table output formats.
    """

    NATIVE = "native"
    CSV = "csv"
    TSV = "tsv"
    XVG = "xvg"



@dataclass(frozen=True)
class AnalysisTable:

    """
    A validated numerical table with scientific column metadata.
    """

    schema: AnalysisSchema
    data: np.ndarray

    def __post_init__(self) -> None:
        data = np.asarray(self.data)

        if data.ndim != 2:
            raise AnalysisOutputError(
                f"Analysis table data must be two-dimensional, got {data.ndim}."
            )

        if data.shape[1] != len(self.schema.columns):
            raise AnalysisOutputError(
                f"Analysis table has {data.shape[1]} data columns but "
                f"{len(self.schema.columns)} metadata columns."
            )

        if data.shape[0] == 0:
            raise AnalysisOutputError("Analysis table contains no data rows.")

        if not np.issubdtype(data.dtype, np.number):
            raise AnalysisOutputError(
                "Analysis table contains non-numeric data."
            )

        object.__setattr__(self, "data", data)

    @classmethod
    def from_columns(
        cls,
        schema: AnalysisSchema,
        columns: Sequence[Sequence[Real] | np.ndarray],
    ) -> "AnalysisTable":
        """
        Build a table from equally sized one-dimensional columns.
        """
        arrays = tuple(np.asarray(column) for column in columns)

        if len(arrays) != len(schema.columns):
            raise AnalysisOutputError(
                f"Analysis table received {len(arrays)} data columns but "
                f"{len(schema.columns)} metadata columns."
            )

        if any(array.ndim != 1 for array in arrays):
            raise AnalysisOutputError(
                "Analysis table columns must be one-dimensional."
            )

        lengths = {len(array) for array in arrays}
        if len(lengths) != 1:
            raise AnalysisOutputError(
                "Analysis table columns must contain the same number of rows."
            )

        return cls(schema=schema, data=np.column_stack(arrays))

    def column(self, field: str) -> np.ndarray:
        """
        Return one data column selected by its stable field identifier.
        """
        try:
            index = self.schema.fields.index(field)
        except ValueError as exception:
            fields = ", ".join(self.schema.fields)
            raise AnalysisOutputError(
                f"Unknown analysis field '{field}'. Available fields: {fields}."
            ) from exception

        return self.data[:, index]



def infer_output_format(filename: str | None) -> AnalysisOutputFormat:
    """
    Infer a requested output format from a recognized filename suffix.

    Unrecognized suffixes and stdout retain the native PQAnalysis format.

    Parameters
    ----------
    filename : str | None
        Requested output filename, or None for stdout.

    Returns
    -------
    AnalysisOutputFormat
        Format selected from the case-insensitive suffix.
    """
    if filename is None:
        return AnalysisOutputFormat.NATIVE

    suffix = Path(filename).suffix.lower()
    return {
        ".csv": AnalysisOutputFormat.CSV,
        ".tsv": AnalysisOutputFormat.TSV,
        ".xvg": AnalysisOutputFormat.XVG,
    }.get(suffix, AnalysisOutputFormat.NATIVE)



def read_analysis_table(filename: str) -> AnalysisTable:
    """
    Read native, CSV, TSV or XVG analysis data by content.

    Parameters
    ----------
    filename : str
        Input table. Its extension is not used for format detection.

    Returns
    -------
    AnalysisTable
        Numeric data and reconstructed scientific schema.

    Raises
    ------
    AnalysisOutputError
        If metadata, headers or numeric rows are malformed.
    """
    path = Path(filename)
    text = path.read_text(encoding="utf-8-sig")
    first_line = next((line for line in text.splitlines() if line.strip()), "")

    if not first_line:
        raise AnalysisOutputError(
            f"Analysis output file '{filename}' is empty."
        )

    if is_xvg(text):
        return _read_xvg_table(text, path)
    if first_line.lstrip().startswith("#"):
        return _read_native_table(text, path)
    if "\t" in first_line:
        return _read_delimited_table(text, path, "\t")
    if "," in first_line:
        return _read_delimited_table(text, path, ",")
    if _is_single_column_delimited(text):
        return _read_delimited_table(text, path, ",")

    return _read_native_table(text, path)



def convert_analysis_output(
    input_file: str,
    output_files: Sequence[str],
    x_field: str | None = None,
    y_fields: Sequence[str] | None = None,
    mode: str | FileWritingMode = "w",
) -> None:
    """
    Convert one analysis table to one or more output files.

    Parameters
    ----------
    input_file : str
        Native, CSV, TSV or PQAnalysis XVG input table.
    output_files : Sequence[str]
        Distinct output paths. Each suffix selects its format.
    x_field : str | None, optional
        XVG x-axis field override, by default None.
    y_fields : Sequence[str] | None, optional
        XVG y-axis fields, by default None.
    mode : str | FileWritingMode, optional
        Output writing mode, by default ``w``.
    """
    if not output_files:
        raise AnalysisOutputError("At least one output file is required.")

    _validate_distinct_paths(input_file, output_files)
    table = read_analysis_table(input_file)
    writers = tuple(
        AnalysisTableWriter(filename, mode=mode) for filename in output_files
    )
    rendered_outputs = tuple(
        (
            writer,
            _render_analysis_table(
                table,
                writer.output_format,
                x_field=x_field,
                y_fields=y_fields,
            ),
        ) for writer in writers
    )

    for writer, rendered in rendered_outputs:
        _write_rendered_output(writer, rendered)



class AnalysisTableWriter(BaseWriter):

    """
    Write an :class:`AnalysisTable` using a suffix-selected format.
    """

    def __init__(
        self,
        filename: str,
        mode: str | FileWritingMode = "w",
    ) -> None:
        """
        Parameters
        ----------
        filename : str
            Output path whose suffix selects the format.
        mode : str | FileWritingMode, optional
            Output writing mode, by default ``w``. Append mode is only
            valid for native text.
        """
        self.output_format = infer_output_format(filename)
        super().__init__(filename, mode=mode)

        if (
            self.output_format != AnalysisOutputFormat.NATIVE and
            self.original_mode == FileWritingMode.APPEND
        ):
            raise AnalysisOutputError(
                f"Appending is not supported for {self.output_format.value} "
                "analysis output."
            )

    def write(
        self,
        table: AnalysisTable,
        x_field: str | None = None,
        y_fields: Sequence[str] | None = None,
    ) -> None:
        """
        Write a table and close the output stream.
        """
        rendered = _render_analysis_table(
            table,
            self.output_format,
            x_field=x_field,
            y_fields=y_fields,
        )
        _write_rendered_output(self, rendered)



class AnalysisDataWriter(BaseWriter):

    """
    Base writer that preserves native formatting and adds table exports.
    """

    schema: AnalysisSchema

    def __init__(
        self,
        filename: str | None,
        mode: str | FileWritingMode = "w",
        export_files: Sequence[str] | None = None,
    ) -> None:
        """
        Parameters
        ----------
        filename : str | None
            Primary output path, or None for native text on stdout.
        mode : str | FileWritingMode, optional
            Output writing mode, by default ``w``.
        export_files : Sequence[str] | None, optional
            Additional distinct output paths, by default None.
        """
        self.output_format = infer_output_format(filename)
        self.export_files = tuple(export_files or ())
        _validate_distinct_paths(filename, self.export_files)
        super().__init__(filename, mode=mode)

        if (
            self.output_format != AnalysisOutputFormat.NATIVE and
            self.original_mode == FileWritingMode.APPEND
        ):
            raise AnalysisOutputError(
                f"Appending is not supported for {self.output_format.value} "
                "analysis output."
            )

        self.export_writers = tuple(
            AnalysisTableWriter(export_file, mode=mode)
            for export_file in self.export_files
        )

    def write_table(
        self,
        table: AnalysisTable,
        native_writer: Callable[[AnalysisTextStream], None],
    ) -> None:
        """
        Write the primary file and all requested additional exports.
        """
        rendered = None
        if self.output_format != AnalysisOutputFormat.NATIVE:
            rendered = _render_analysis_table(table, self.output_format)

        super().open()
        try:
            if rendered is None:
                native_writer(self.file)
            else:
                self.file.write(rendered)
        finally:
            super().close()

        for writer in self.export_writers:
            writer.write(table)



def write_analysis_table(
    table: AnalysisTable,
    file: AnalysisTextStream,
    output_format: AnalysisOutputFormat,
    x_field: str | None = None,
    y_fields: Sequence[str] | None = None,
) -> None:
    """
    Write an analysis table to an open text stream.

    Parameters
    ----------
    table : AnalysisTable
        Validated data and scientific schema.
    file : TextIO
        Open UTF-8 text stream.
    output_format : AnalysisOutputFormat
        Format to write.
    x_field : str | None, optional
        XVG x-axis field override, by default None.
    y_fields : Sequence[str] | None, optional
        XVG y-axis fields, by default None.
    """
    if output_format == AnalysisOutputFormat.NATIVE:
        _write_native_table(table, file)
    elif output_format == AnalysisOutputFormat.CSV:
        _write_delimited_table(table, file, ",")
    elif output_format == AnalysisOutputFormat.TSV:
        _write_delimited_table(table, file, "\t")
    elif output_format == AnalysisOutputFormat.XVG:
        _write_xvg_table(table, file, x_field=x_field, y_fields=y_fields)
    else:
        raise AnalysisOutputError(
            f"Unsupported analysis output format '{output_format}'."
        )



def _render_analysis_table(
    table: AnalysisTable,
    output_format: AnalysisOutputFormat,
    x_field: str | None = None,
    y_fields: Sequence[str] | None = None,
) -> str:
    buffer = io.StringIO()
    write_analysis_table(
        table,
        buffer,
        output_format,
        x_field=x_field,
        y_fields=y_fields,
    )
    return buffer.getvalue()



def _write_rendered_output(
    writer: AnalysisTableWriter,
    rendered: str,
) -> None:
    writer.open()
    try:
        writer.file.write(rendered)
    finally:
        writer.close()



def _is_single_column_delimited(text: str) -> bool:
    lines = tuple(line.strip() for line in text.splitlines() if line.strip())
    if len(lines) < 2 or lines[0].startswith("#"):
        return False

    try:
        float(lines[0])
        return False
    except ValueError:
        pass

    for line in lines[1:]:
        try:
            float(line)
        except ValueError:
            return False
    return True



def _read_xvg_table(text: str, path: Path) -> AnalysisTable:
    try:
        document = read_xvg(text, path)
    except XVGFormatError as exception:
        raise AnalysisOutputError(str(exception)) from exception

    schema = _schema_from_metadata(
        title=document.title,
        fields=document.fields,
        symbols=document.symbols,
        units=document.units,
    )
    schema = replace(
        schema,
        plot=_plot_preset_for_fields(
            schema,
            document.x_field,
            document.y_fields,
        ),
    )
    return AnalysisTable(schema=schema, data=document.data)



def _read_native_table(text: str, path: Path) -> AnalysisTable:
    metadata, numeric_lines = _scan_native_lines(text)

    data = _load_native_data(numeric_lines, path)
    fields, symbols, units = _native_column_metadata(metadata, data.shape[1])

    _validate_metadata(fields, symbols, units, data.shape[1], path)
    schema = _schema_from_metadata(
        title=metadata.get("title", path.stem),
        fields=fields,
        symbols=symbols,
        units=units,
    )
    return AnalysisTable(schema=schema, data=data)



def _scan_native_lines(text: str) -> tuple[dict[str, str], list[str]]:
    metadata: dict[str, str] = {}
    numeric_lines = []

    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if not stripped.startswith("#"):
            numeric_lines.append(line)
            continue

        content = stripped[1:].strip()
        if content.startswith("PQAnalysis:"):
            metadata["title"] = content.removeprefix("PQAnalysis:").strip()
            continue

        key, _, value = content.partition(" ")
        if key in {"FIELDS", "SYMBOLS", "UNITS"}:
            metadata[key.lower()] = value.strip()

    return metadata, numeric_lines



def _load_native_data(numeric_lines: Sequence[str], path: Path) -> np.ndarray:
    if not numeric_lines:
        raise AnalysisOutputError(
            f"Analysis output file '{path}' contains no numeric rows."
        )

    try:
        return np.loadtxt(io.StringIO("\n".join(numeric_lines)), ndmin=2)
    except ValueError as exception:
        raise AnalysisOutputError(
            f"Analysis output file '{path}' contains invalid numeric data."
        ) from exception



def _native_column_metadata(
    metadata: dict[str, str],
    n_columns: int,
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    fields = tuple(metadata.get("fields", "").split())
    symbols = tuple(metadata.get("symbols", "").split())
    units = tuple(metadata.get("units", "").split())

    if not fields:
        fields = tuple(f"column_{index}" for index in range(1, n_columns + 1))

    return fields, symbols or fields, units or tuple("1" for _ in fields)



def _read_delimited_table(
    text: str, path: Path, delimiter: str
) -> AnalysisTable:
    rows = [
        row for row in csv.reader(io.StringIO(text), delimiter=delimiter)
        if row and any(value.strip() for value in row)
    ]

    if len(rows) < 2:
        raise AnalysisOutputError(
            f"Delimited analysis output file '{path}' contains no data rows."
        )

    fields = tuple(value.strip() for value in rows[0])
    _validate_fields(fields, path)

    if any(len(row) != len(fields) for row in rows[1:]):
        raise AnalysisOutputError(
            f"Delimited analysis output file '{path}' has inconsistent row widths."
        )

    try:
        data = np.asarray(
            [[float(value) for value in row] for row in rows[1:]],
            dtype=float,
        )
    except ValueError as exception:
        raise AnalysisOutputError(
            f"Delimited analysis output file '{path}' contains invalid numeric data."
        ) from exception

    schema = _schema_for_fields(fields, path.stem)
    return AnalysisTable(schema=schema, data=data)



def _schema_from_metadata(
    title: str,
    fields: tuple[str, ...],
    symbols: tuple[str, ...],
    units: tuple[str, ...],
) -> AnalysisSchema:
    known = _schema_for_fields(fields, title)
    columns = tuple(
        AnalysisColumn(field, symbol, unit)
        for field, symbol, unit in zip(fields, symbols, units)
    )
    return replace(known, title=title, columns=columns)



def _schema_for_fields(fields: tuple[str, ...], title: str) -> AnalysisSchema:
    if fields in KNOWN_SCHEMAS:
        return KNOWN_SCHEMAS[fields]

    expected_mode_fields = tuple(
        f"mode_{index}" for index in range(1, len(fields) + 1)
    )
    if fields == expected_mode_fields:
        return normal_mode_schema(len(fields))

    columns = tuple(AnalysisColumn(field, field, "1") for field in fields)
    if len(fields) == 1:
        plot = AnalysisPlotPreset(
            x_field=None,
            y_fields=(fields[0], ),
            x_label="Row",
            y_label=fields[0],
            legends=(fields[0], ),
        )
    else:
        plot = AnalysisPlotPreset(
            x_field=fields[0],
            y_fields=(fields[1], ),
            x_label=fields[0],
            y_label=fields[1],
            legends=(fields[1], ),
        )

    return AnalysisSchema(title=title, columns=columns, plot=plot)



def _validate_metadata(
    fields: tuple[str, ...],
    symbols: tuple[str, ...],
    units: tuple[str, ...],
    n_columns: int,
    path: Path,
) -> None:
    counts = (len(fields), len(symbols), len(units), n_columns)
    if len(set(counts)) != 1:
        raise AnalysisOutputError(
            f"Analysis output file '{path}' has inconsistent FIELDS, SYMBOLS, "
            "UNITS and numeric column counts."
        )
    _validate_fields(fields, path)



def _validate_fields(fields: tuple[str, ...], path: Path) -> None:
    if not fields or any(not field for field in fields):
        raise AnalysisOutputError(
            f"Analysis output file '{path}' contains a blank field name."
        )
    if len(set(fields)) != len(fields):
        raise AnalysisOutputError(
            f"Analysis output file '{path}' contains duplicate field names."
        )



def _write_native_table(
    table: AnalysisTable, file: AnalysisTextStream
) -> None:
    print(f"# PQAnalysis: {table.schema.title}", file=file)
    print(f"# FIELDS {' '.join(table.schema.fields)}", file=file)
    print(
        f"# SYMBOLS {' '.join(column.symbol for column in table.schema.columns)}",
        file=file,
    )
    print(
        f"# UNITS {' '.join(column.unit for column in table.schema.columns)}",
        file=file,
    )

    for row in table.data:
        print(" ".join(_format_number(value) for value in row), file=file)



def _write_delimited_table(
    table: AnalysisTable,
    file: AnalysisTextStream,
    delimiter: str,
) -> None:
    writer = csv.writer(file, delimiter=delimiter, lineterminator="\n")
    writer.writerow(table.schema.fields)
    writer.writerows(
        [_format_number(value) for value in row] for row in table.data
    )



def _write_xvg_table(
    table: AnalysisTable,
    file: AnalysisTextStream,
    x_field: str | None = None,
    y_fields: Sequence[str] | None = None,
) -> None:
    preset = table.schema.plot

    if preset is None:
        raise AnalysisOutputError(
            f"Analysis table '{table.schema.title}' has no xmgrace plot preset."
        )

    selected_x = x_field if x_field is not None else preset.x_field
    selected_y = tuple(y_fields) if y_fields is not None else preset.y_fields
    _validate_xvg_projection(table.schema, selected_x, selected_y)

    plot = _plot_preset_for_fields(table.schema, selected_x, selected_y)
    x_values = (
        np.arange(1, table.data.shape[0] +
                  1) if selected_x is None else table.column(selected_x)
    )

    selected_legends = dict(zip(plot.y_fields, plot.legends))
    legends = tuple(
        selected_legends.get(field, _column_label(table.schema, field))
        for field in table.schema.fields
    )
    write_xvg(
        file=file,
        document=XVGDocument(
            title=table.schema.title,
            fields=table.schema.fields,
            symbols=tuple(column.symbol for column in table.schema.columns),
            units=tuple(column.unit for column in table.schema.columns),
            x_field=selected_x,
            y_fields=selected_y,
            data=table.data,
        ),
        style=XVGPlotStyle(
            x_label=plot.x_label,
            y_label=plot.y_label,
            legends=legends,
            x_values=x_values,
        ),
    )



def _validate_xvg_projection(
    schema: AnalysisSchema,
    x_field: str | None,
    y_fields: tuple[str, ...],
) -> None:
    if x_field is not None:
        _column_label(schema, x_field)
    if not y_fields:
        raise AnalysisOutputError(
            "At least one y field is required for XVG output."
        )
    if len(set(y_fields)) != len(y_fields):
        raise AnalysisOutputError("XVG y fields must be distinct.")
    for field in y_fields:
        _column_label(schema, field)



def _plot_preset_for_fields(
    schema: AnalysisSchema,
    x_field: str | None,
    y_fields: tuple[str, ...],
) -> AnalysisPlotPreset:
    preset = schema.plot
    if preset is not None and (
        x_field == preset.x_field and y_fields == preset.y_fields
    ):
        return preset

    if x_field is None:
        x_label = (
            preset.x_label
            if preset is not None and preset.x_field is None else "Row"
        )
    else:
        x_label = _column_label(schema, x_field)

    y_label = (
        _column_label(schema, y_fields[0]) if len(y_fields) == 1 else "Value"
    )
    legends = tuple(_column_label(schema, field) for field in y_fields)
    return AnalysisPlotPreset(
        x_field=x_field,
        y_fields=y_fields,
        x_label=x_label,
        y_label=y_label,
        legends=legends,
    )



def _column_label(schema: AnalysisSchema, field: str) -> str:
    try:
        column = next(
            column for column in schema.columns if column.field == field
        )
    except StopIteration as exception:
        fields = ", ".join(schema.fields)
        raise AnalysisOutputError(
            f"Unknown analysis field '{field}'. Available fields: {fields}."
        ) from exception

    symbol = _unicode_to_grace(column.symbol)
    if column.unit == "1":
        return symbol
    return f"{symbol} [{_unicode_to_grace(column.unit)}]"



def _unicode_to_grace(value: str) -> str:
    replacements = {
        "Å": r"\cE\C",
        "²": r"\S2\N",
        "³": r"\S3\N",
        "⁻¹": r"\S-1\N",
        "−": "-",
        "∕": "/",
        "·": " ",
        "‖": "|",
        "⟨": "<",
        "⟩": ">",
        "Δ": "Delta",
        "μ": "mu",
        "ν̃": "nu~",
        "Ĉ": "C^",
        "ᵢ": r"\si\N",
        "ⱼ": r"\sj\N",
        "ᴵᴿ": r"\SIR\N",
        "α": "alpha",
    }
    converted = value
    for source, target in replacements.items():
        converted = converted.replace(source, target)
    return converted



def _format_number(value: object) -> str:
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return repr(float(value))



def _validate_distinct_paths(
    primary_file: str | None,
    output_files: Sequence[str],
) -> None:
    paths = []
    if primary_file is not None:
        paths.append(Path(primary_file).expanduser().resolve())
    paths.extend(Path(path).expanduser().resolve() for path in output_files)

    if len(paths) != len(set(paths)):
        raise AnalysisOutputError("Analysis output paths must be distinct.")
