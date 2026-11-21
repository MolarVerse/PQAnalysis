"""
Read and write reversible PQAnalysis XVG documents.
"""

import io
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

import numpy as np

XVG_VERSION = "1"
XVG_ROW_INDEX = "@row"
XVGTextStream = TextIO | io.StringIO
_TARGET_PATTERN = re.compile(r"target\s+G(\d+)\.S(\d+)", re.IGNORECASE)
# Grace colours: blue, red, green4, orange, violet, maroon, turquoise,
# indigo, magenta, brown, cyan, black (yellow and grey skipped).
_GRACE_COLORS = (4, 2, 15, 11, 8, 13, 14, 12, 10, 6, 9, 1)
_GRACE_LINE_STYLES = (1, 3, 2, 4)
_VIEW_LEFT, _VIEW_RIGHT = 0.2, 0.95
_VIEW_BOTTOM, _VIEW_TOP = 0.12, 0.88
_GRAPH_GAP = 0.06



class XVGFormatError(ValueError):

    """Raised when an XVG document cannot reconstruct an analysis table."""



@dataclass(frozen=True)
class XVGDocument:  # pylint: disable=too-many-instance-attributes

    """
    Complete analysis-table data and plot fields decoded from XVG.

    Non-empty ``graphs`` give the y fields of each stacked Grace graph.
    """

    title: str
    fields: tuple[str, ...]
    symbols: tuple[str, ...]
    units: tuple[str, ...]
    x_field: str | None
    y_fields: tuple[str, ...]
    data: np.ndarray
    graphs: tuple[tuple[str, ...], ...] = ()



@dataclass(frozen=True)
class XVGPlotStyle:

    """
    Grace labels, legends and x-axis values for an XVG document.

    ``graph_labels`` are the y labels of the document graphs.
    """

    x_label: str
    y_label: str
    legends: tuple[str, ...]
    x_values: np.ndarray
    graph_labels: tuple[str, ...] = ()



@dataclass(frozen=True)
class _XVGScan:

    metadata: dict[str, str]
    targets: tuple[tuple[int, int], ...]
    graphs: tuple[tuple[str, ...], ...]
    set_types: tuple[str | None, ...]
    numeric_sets: tuple[tuple[str, ...], ...]
    malformed_layout: bool



@dataclass
class _XVGScanState:

    metadata: dict[str, str] = field(default_factory=dict)
    targets: list[tuple[int, int]] = field(default_factory=list)
    graphs: list[tuple[str, ...]] = field(default_factory=list)
    set_types: list[str | None] = field(default_factory=list)
    numeric_sets: list[list[str]] = field(default_factory=list)
    current_set: int | None = None
    malformed_layout: bool = False



def is_xvg(text: str) -> bool:
    """Return whether text contains Grace directives or a PQAnalysis marker."""
    return any(
        line.lstrip().startswith(("@", "# PQAnalysis-XVG:"))
        for line in text.splitlines()
    )



def read_xvg(text: str, path: Path) -> XVGDocument:
    """Decode a reversible PQAnalysis Grace XY-block document."""
    scan = _scan_xvg_lines(text)

    if scan.metadata.get("xvg_version") != XVG_VERSION:
        raise XVGFormatError(
            f"XVG file '{path}' does not contain reversible PQAnalysis "
            "table metadata."
        )
    if scan.malformed_layout:
        raise XVGFormatError(
            f"PQAnalysis XVG file '{path}' must use complete XY data-set "
            "blocks."
        )

    fields = tuple(scan.metadata.get("fields", "").split())
    symbols = tuple(scan.metadata.get("symbols", "").split())
    units = tuple(scan.metadata.get("units", "").split())
    _validate_columns(fields, symbols, units, len(scan.numeric_sets), path)
    _validate_set_layout(scan, len(fields), path)
    axis, table_data = _load_xvg_sets(scan.numeric_sets, path)

    x_field = _read_x_field(scan.metadata, fields, path)
    y_fields = tuple(scan.metadata.get("xvg_y_fields", "").split())
    _validate_plot_fields(fields, y_fields, path)
    _validate_graphs(scan.graphs, y_fields, path)
    _validate_axis(axis, table_data, fields, x_field, path)

    return XVGDocument(
        title=scan.metadata.get("title", path.stem),
        fields=fields,
        symbols=symbols,
        units=units,
        x_field=x_field,
        y_fields=y_fields,
        data=table_data,
        graphs=scan.graphs,
    )



def write_xvg(
    file: XVGTextStream,
    document: XVGDocument,
    style: XVGPlotStyle,
) -> None:
    """Write a complete analysis table as standard Grace XY data sets."""
    _validate_write_data(document, style)
    _write_metadata(file, document)
    if document.graphs:
        _write_graphs(file, document, style)
        return

    print(f'@    title "{_quote(document.title)}"', file=file)
    print(f'@    xaxis label "{_quote(style.x_label)}"', file=file)
    print(f'@    yaxis label "{_quote(style.y_label)}"', file=file)
    print("@    legend on", file=file)

    columns = document.data.T
    for index, (field_name, legend, values) in enumerate(
        zip(document.fields, style.legends, columns)
    ):
        _write_set_data(file, (0, index), style.x_values, values)
        print(f'@    s{index} legend "{_quote(legend)}"', file=file)
        if field_name not in document.y_fields:
            print(f"@    s{index} hidden true", file=file)



def _write_metadata(file: XVGTextStream, document: XVGDocument) -> None:
    print(f"# PQAnalysis-XVG: {XVG_VERSION}", file=file)
    print("# Generated by PQAnalysis", file=file)
    print(f"# PQAnalysis: {document.title}", file=file)
    print(f"# FIELDS {' '.join(document.fields)}", file=file)
    print(f"# SYMBOLS {' '.join(document.symbols)}", file=file)
    print(f"# UNITS {' '.join(document.units)}", file=file)
    print(
        f"# XVG_X_FIELD {document.x_field or XVG_ROW_INDEX}",
        file=file,
    )
    print(f"# XVG_Y_FIELDS {' '.join(document.y_fields)}", file=file)
    for graph in document.graphs:
        print(f"# XVG_GRAPH {' '.join(graph)}", file=file)



def _write_set_data(
    file: XVGTextStream,
    target: tuple[int, int],
    x_values: np.ndarray,
    values: np.ndarray,
) -> None:
    print(f"@target G{target[0]}.S{target[1]}", file=file)
    print("@type xy", file=file)
    for x_value, value in zip(x_values, values):
        print(f"{_format_number(x_value)} {_format_number(value)}", file=file)
    print("&", file=file)



def _graph_targets(
    fields: tuple[str, ...],
    graphs: tuple[tuple[str, ...], ...],
) -> tuple[tuple[int, int], ...]:
    graph_of = {
        name: index
        for index, graph in enumerate(graphs)
        for name in graph
    }
    counts = [0] * len(graphs)
    targets = []
    for name in fields:
        graph = graph_of.get(name, 0)
        targets.append((graph, counts[graph]))
        counts[graph] += 1
    return tuple(targets)



def _write_graphs(
    file: XVGTextStream,
    document: XVGDocument,
    style: XVGPlotStyle,
) -> None:
    targets = _graph_targets(document.fields, document.graphs)
    plotted = tuple(name for graph in document.graphs for name in graph)
    style_order = plotted + tuple(
        name for name in document.fields if name not in plotted
    )
    print("@version 50125", file=file)
    for index in range(len(document.graphs)):
        print(f"@g{index} on", file=file)
        print(f"@g{index} hidden false", file=file)
        print(f"@g{index} type XY", file=file)

    for index, graph in enumerate(document.graphs):
        values = document.data[:,
                               [document.fields.index(name) for name in graph]]
        _write_graph_header(file, index, document, style, values)
        for target, name, legend in zip(
            targets, document.fields, style.legends
        ):
            if target[0] != index:
                continue
            _write_set_style(file, target[1], style_order.index(name), legend)
            if name not in graph:
                print(f"@    s{target[1]} hidden true", file=file)

    for target, values in zip(targets, document.data.T):
        _write_set_data(file, target, style.x_values, values)
    for index in range(len(document.graphs)):
        print(f"@with g{index}", file=file)
        print("@autoscale", file=file)



def _write_graph_header(
    file: XVGTextStream,
    index: int,
    document: XVGDocument,
    style: XVGPlotStyle,
    values: np.ndarray,
) -> None:
    n_graphs = len(document.graphs)
    height = (
        _VIEW_TOP - _VIEW_BOTTOM - (n_graphs - 1) * _GRAPH_GAP
    ) / n_graphs
    top = _VIEW_TOP - index * (height + _GRAPH_GAP)
    label_size = 1.2 if n_graphs == 1 else 1.0
    tick_size = 1.0 if n_graphs == 1 else 0.85

    print(f"@with g{index}", file=file)
    print(
        f"@    view {_VIEW_LEFT:.4f}, {top - height:.4f}, "
        f"{_VIEW_RIGHT:.4f}, {top:.4f}",
        file=file,
    )
    if index == 0:
        print(f'@    title "{_quote(document.title)}"', file=file)
        print("@    title size 1.2", file=file)
    if index == n_graphs - 1:
        print(f'@    xaxis label "{_quote(style.x_label)}"', file=file)
        print(f"@    xaxis label char size {label_size}", file=file)
        print(f"@    xaxis ticklabel char size {tick_size}", file=file)
    else:
        print("@    xaxis ticklabel off", file=file)
    y_label = style.graph_labels[index]
    print(f'@    yaxis label "{_quote(y_label)}"', file=file)
    print(f"@    yaxis label char size {label_size}", file=file)
    print(f"@    yaxis ticklabel char size {tick_size}", file=file)
    precision = _tick_precision(values)
    if precision is not None:
        print("@    yaxis ticklabel format decimal", file=file)
        print(f"@    yaxis ticklabel prec {precision}", file=file)
    print("@    legend on", file=file)
    print("@    legend loctype view", file=file)
    print(f"@    legend {_VIEW_RIGHT + 0.02:.4f}, {top:.4f}", file=file)
    print("@    legend box linestyle 0", file=file)
    print("@    legend char size 0.8", file=file)



def _tick_precision(values: np.ndarray) -> int | None:
    """Return decimal tick digits for the data range, None for Grace's."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return None
    span = float(np.ptp(finite))
    magnitude = float(np.max(np.abs(finite)))
    if span <= 0 or magnitude >= 1e9:
        return None
    precision = max(0, int(np.ceil(-np.log10(span / 5))))
    return precision if precision <= 6 else None



def _write_set_style(
    file: XVGTextStream,
    set_index: int,
    style_index: int,
    legend: str,
) -> None:
    color = _GRACE_COLORS[style_index % len(_GRACE_COLORS)]
    line_style = _GRACE_LINE_STYLES[(style_index // len(_GRACE_COLORS)) %
                                    len(_GRACE_LINE_STYLES)]
    print(f"@    s{set_index} line type 1", file=file)
    print(f"@    s{set_index} line linestyle {line_style}", file=file)
    print(f"@    s{set_index} line linewidth 2.0", file=file)
    print(f"@    s{set_index} line color {color}", file=file)
    print(f'@    s{set_index} legend "{_quote(legend)}"', file=file)



def _validate_graphs(
    graphs: tuple[tuple[str, ...], ...],
    y_fields: tuple[str, ...],
    path: Path,
) -> None:
    graph_fields = tuple(name for graph in graphs for name in graph)
    if graphs and (
        not all(graphs) or sorted(graph_fields) != sorted(y_fields)
    ):
        raise XVGFormatError(
            f"PQAnalysis XVG file '{path}' has graphs that do not match "
            "its y fields."
        )



def _validate_write_data(document: XVGDocument, style: XVGPlotStyle) -> None:
    column_counts = (
        len(document.fields),
        len(document.symbols),
        len(document.units),
        len(style.legends),
    )
    if document.data.ndim != 2 or len(set(column_counts)) != 1:
        raise XVGFormatError("XVG table metadata and data columns must match.")
    if document.data.shape[1] != len(document.fields):
        raise XVGFormatError("XVG table metadata and data columns must match.")
    if document.data.shape[0] != len(style.x_values):
        raise XVGFormatError("XVG x-axis and table row counts must match.")
    if len(style.graph_labels) != len(document.graphs):
        raise XVGFormatError("XVG graphs and graph labels must match.")



def _scan_xvg_lines(text: str, ) -> _XVGScan:
    state = _XVGScanState()

    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        _scan_xvg_line(stripped, line, state)

    if state.current_set is not None:
        state.malformed_layout = True

    return _XVGScan(
        metadata=state.metadata,
        targets=tuple(state.targets),
        graphs=tuple(state.graphs),
        set_types=tuple(state.set_types),
        numeric_sets=tuple(tuple(lines) for lines in state.numeric_sets),
        malformed_layout=state.malformed_layout,
    )



def _scan_xvg_line(stripped: str, line: str, state: _XVGScanState) -> None:
    if stripped.startswith("# XVG_GRAPH "):
        state.graphs.append(tuple(stripped.split()[2:]))
    elif stripped.startswith("#"):
        _read_metadata_line(stripped[1:].strip(), state.metadata)
    elif stripped == "&":
        if state.current_set is None:
            state.malformed_layout = True
        state.current_set = None
    elif stripped.startswith("@"):
        _read_xvg_directive(stripped[1:].strip(), state)
    elif state.current_set is None:
        state.malformed_layout = True
    else:
        state.numeric_sets[state.current_set].append(line)



def _read_xvg_directive(directive: str, state: _XVGScanState) -> None:
    target_match = _TARGET_PATTERN.fullmatch(directive)
    if target_match is not None:
        if state.current_set is not None:
            state.malformed_layout = True
        state.targets.append(
            (int(target_match.group(1)), int(target_match.group(2)))
        )
        state.set_types.append(None)
        state.numeric_sets.append([])
        state.current_set = len(state.numeric_sets) - 1
        return

    if directive.lower().startswith("target "):
        state.malformed_layout = True
    elif directive.lower().startswith("type "):
        if state.current_set is None:
            state.malformed_layout = True
        else:
            state.set_types[state.current_set] = directive.split(None, 1)[1]



def _read_metadata_line(content: str, metadata: dict[str, str]) -> None:
    if content.startswith("PQAnalysis-XVG:"):
        metadata["xvg_version"] = content.removeprefix("PQAnalysis-XVG:"
                                                       ).strip()
        return
    if content.startswith("PQAnalysis:"):
        metadata["title"] = content.removeprefix("PQAnalysis:").strip()
        return

    key, _, value = content.partition(" ")
    metadata_keys = {
        "FIELDS": "fields",
        "SYMBOLS": "symbols",
        "UNITS": "units",
        "XVG_X_FIELD": "xvg_x_field",
        "XVG_Y_FIELDS": "xvg_y_fields",
    }
    if key in metadata_keys:
        metadata[metadata_keys[key]] = value.strip()



def _validate_set_layout(scan: _XVGScan, n_columns: int, path: Path) -> None:
    next_set: dict[int, int] = {}
    sequential = len(scan.targets) == n_columns
    for graph, set_index in scan.targets:
        sequential = sequential and set_index == next_set.get(graph, 0)
        next_set[graph] = set_index + 1
    if not sequential or (0 not in next_set and n_columns):
        raise XVGFormatError(
            f"PQAnalysis XVG file '{path}' must define sequential Grace sets."
        )
    if any(
        set_type is None or set_type.lower() != "xy"
        for set_type in scan.set_types
    ):
        raise XVGFormatError(
            f"PQAnalysis XVG file '{path}' must use XY data sets."
        )



def _load_xvg_sets(
    numeric_sets: tuple[tuple[str, ...], ...],
    path: Path,
) -> tuple[np.ndarray, np.ndarray]:
    loaded_sets = []

    for numeric_lines in numeric_sets:
        if not numeric_lines:
            raise XVGFormatError(
                f"PQAnalysis XVG file '{path}' contains an empty data set."
            )
        try:
            data_set = np.loadtxt(numeric_lines, ndmin=2)
        except ValueError as exception:
            raise XVGFormatError(
                f"PQAnalysis XVG file '{path}' contains invalid numeric data."
            ) from exception
        if data_set.shape[1] != 2:
            raise XVGFormatError(
                f"PQAnalysis XVG file '{path}' data sets must contain x and y."
            )
        loaded_sets.append(data_set)

    axis = loaded_sets[0][:, 0]
    if any(
        data_set.shape[0] != axis.shape[0] or
        not np.allclose(data_set[:, 0], axis, equal_nan=True)
        for data_set in loaded_sets[1:]
    ):
        raise XVGFormatError(
            f"PQAnalysis XVG file '{path}' has inconsistent x-axis columns."
        )
    table_data = np.column_stack([data_set[:, 1] for data_set in loaded_sets])
    return axis, table_data



def _validate_columns(
    fields: tuple[str, ...],
    symbols: tuple[str, ...],
    units: tuple[str, ...],
    n_columns: int,
    path: Path,
) -> None:
    counts = (len(fields), len(symbols), len(units), n_columns)
    if len(set(counts)) != 1:
        raise XVGFormatError(
            f"Analysis output file '{path}' has inconsistent FIELDS, SYMBOLS, "
            "UNITS and numeric column counts."
        )
    if not fields or any(not field for field in fields):
        raise XVGFormatError(
            f"Analysis output file '{path}' contains a blank field name."
        )
    if len(set(fields)) != len(fields):
        raise XVGFormatError(
            f"Analysis output file '{path}' contains duplicate field names."
        )



def _read_x_field(
    metadata: dict[str, str],
    fields: tuple[str, ...],
    path: Path,
) -> str | None:
    value = metadata.get("xvg_x_field")
    if value is None:
        raise XVGFormatError(
            f"PQAnalysis XVG file '{path}' does not define XVG_X_FIELD."
        )
    if value == XVG_ROW_INDEX:
        return None
    if value not in fields:
        raise XVGFormatError(
            f"PQAnalysis XVG file '{path}' uses unknown x field '{value}'."
        )
    return value



def _validate_plot_fields(
    fields: tuple[str, ...],
    y_fields: tuple[str, ...],
    path: Path,
) -> None:
    if not y_fields:
        raise XVGFormatError(
            f"PQAnalysis XVG file '{path}' defines no y fields."
        )
    if len(set(y_fields)) != len(y_fields):
        raise XVGFormatError(
            f"PQAnalysis XVG file '{path}' defines duplicate y fields."
        )

    unknown = tuple(field for field in y_fields if field not in fields)
    if unknown:
        raise XVGFormatError(
            f"PQAnalysis XVG file '{path}' uses unknown y field "
            f"'{unknown[0]}'."
        )



def _validate_axis(
    axis: np.ndarray,
    table_data: np.ndarray,
    fields: tuple[str, ...],
    x_field: str | None,
    path: Path,
) -> None:
    if x_field is None:
        expected = np.arange(1, table_data.shape[0] + 1)
    else:
        expected = table_data[:, fields.index(x_field)]

    if not np.allclose(axis, expected, equal_nan=True):
        raise XVGFormatError(
            f"PQAnalysis XVG file '{path}' has an inconsistent x-axis column."
        )



def _quote(value: str) -> str:
    return value.replace('"', r'\"')



def _format_number(value: object) -> str:
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return repr(float(value))
