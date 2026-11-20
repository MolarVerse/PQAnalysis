"""
Tests for converting PQ and QMCFC energy files to analysis tables.
"""

import csv
import io
import shutil
from pathlib import Path

import numpy as np
import pytest

from PQAnalysis.analysis._energy_table import (
    MAX_DEFAULT_GRAPHS,
    default_energy_plot,
    energy_plot_preset,
    grace_unit,
)
from PQAnalysis.analysis._xvg import (
    XVGDocument,
    XVGFormatError,
    XVGPlotStyle,
    _tick_precision,
    write_xvg,
)
from PQAnalysis.analysis.output import (
    AnalysisColumn,
    AnalysisOutputError,
    AnalysisSchema,
    EnergyFileOptions,
    convert_analysis_output,
    read_analysis_table,
    read_energy_table,
)
from PQAnalysis.traj import MDEngineFormat

DATA_DIR = Path(__file__).parents[1] / "data" / "convertEnergy"

PQ_FIELDS = (
    "SIMULATION-TIME",
    "TEMPERATURE",
    "PRESSURE",
    "E(TOT)",
    "E(QM)",
    "N(QM-ATOMS)",
    "E(KIN)",
    "E(INTRA)",
    "VOLUME",
    "DENSITY",
    "MOMENTUM",
    "LOOPTIME",
)
PQ_CSV_HEADER = (
    "SIMULATION-TIME [step],TEMPERATURE [K],PRESSURE [bar],"
    "E(TOT) [kcal/mol],E(QM) [kcal/mol],N(QM-ATOMS),E(KIN) [kcal/mol],"
    "E(INTRA) [kcal/mol],VOLUME [A^3],DENSITY [g/cm^3],"
    "MOMENTUM [amuA/fs],LOOPTIME [s]"
)



@pytest.fixture
def energy_dir(tmp_path):
    """Copy the energy fixtures to a writable directory."""
    for path in DATA_DIR.iterdir():
        shutil.copy(path, tmp_path / path.name)
    return tmp_path



def test_pq_energy_table_uses_info_names_and_units(energy_dir):
    table = read_energy_table(str(energy_dir / "pq.en"))

    assert table.schema.fields == PQ_FIELDS
    assert table.schema.title == "PQ energy file pq"
    units = {column.field: column.unit for column in table.schema.columns}
    assert units["SIMULATION-TIME"] == "step"
    assert units["TEMPERATURE"] == "K"
    assert units["N(QM-ATOMS)"] == "1"
    assert units["VOLUME"] == "A^3"
    assert table.data.shape == (2, 12)
    np.testing.assert_array_equal(table.column("SIMULATION-TIME"), [1, 2])

    plot = table.schema.plot
    assert plot.x_field == "SIMULATION-TIME"
    assert plot.y_fields == ("E(TOT)", "E(QM)", "E(KIN)")
    assert plot.x_label == "Simulation time [step]"
    assert plot.legends == plot.y_fields
    assert [graph.y_fields for graph in plot.graphs] == [
        ("E(TOT)", ), ("E(QM)", ), ("E(KIN)", )
    ]
    assert plot.graphs[0].y_label == "E(TOT) [kcal/mol]"



def test_pq_energy_timestep_gives_time_in_ps(energy_dir):
    table = read_energy_table(
        str(energy_dir / "pq.en"), EnergyFileOptions(timestep_fs=2.0)
    )

    assert table.schema.columns[0].unit == "ps"
    assert table.schema.plot.x_label == "Simulation time [ps]"
    np.testing.assert_allclose(table.column("SIMULATION-TIME"), [0.002, 0.004])



def test_qmcfc_energy_table(energy_dir):
    table = read_energy_table(
        str(energy_dir / "qmcfc.en"),
        EnergyFileOptions(engine_format=MDEngineFormat.QMCFC),
    )

    assert table.schema.fields[:4] == (
        "SIMULATION_TIME",
        "QM_MOLECULES",
        "TEMPERATURE",
        "PRESSURE",
    )
    assert len(table.schema.fields) == 24
    assert {column.unit for column in table.schema.columns} == {"1"}
    plot = table.schema.plot
    assert plot.x_field == "SIMULATION_TIME"
    assert plot.x_label == "Simulation time"
    assert plot.y_fields == ("E(QM)", "E(MM)", "E(KIN)", "E(INTRA)")
    assert plot.graphs[0].y_label == "E(QM)"



def test_energy_table_without_info(energy_dir):
    (energy_dir / "pq.info").unlink()

    table = read_energy_table(str(energy_dir / "pq.en"))

    assert table.schema.fields[:2] == ("column_1", "column_2")
    plot = table.schema.plot
    assert plot.x_field == "column_1"
    assert len(plot.y_fields) == 11
    assert len(plot.graphs) == 1
    assert plot.graphs[0].y_label == "Value"



def test_energy_table_with_explicit_info_file(energy_dir):
    shutil.copy(energy_dir / "pq.en", energy_dir / "other.dat")

    table = read_energy_table(
        str(energy_dir / "other.dat"),
        EnergyFileOptions(info_file=str(energy_dir / "pq.info")),
    )

    assert table.schema.fields == PQ_FIELDS



def test_convert_energy_to_csv_tsv_xvg_and_native(energy_dir):
    outputs = [
        str(energy_dir / f"pq.{suffix}")
        for suffix in ("csv", "tsv", "xvg", "dat")
    ]

    convert_analysis_output(str(energy_dir / "pq.en"), outputs)

    csv_lines = (energy_dir / "pq.csv").read_text(encoding="utf-8"
                                                  ).splitlines()
    assert csv_lines[0] == PQ_CSV_HEADER
    assert csv_lines[1].startswith("1.0,298.066019722598,")
    tsv_header = (energy_dir / "pq.tsv").read_text(encoding="utf-8")
    assert tsv_header.startswith(PQ_CSV_HEADER.replace(",", "\t") + "\n")
    native = (energy_dir / "pq.dat").read_text(encoding="utf-8")
    assert f"# FIELDS {' '.join(PQ_FIELDS)}\n" in native

    xvg = (energy_dir / "pq.xvg").read_text(encoding="utf-8")
    assert '@    xaxis label "Simulation time [step]"' in xvg
    assert '@    yaxis label "E(TOT) [kcal/mol]"' in xvg
    assert '@    s1 legend "TEMPERATURE [K]"' in xvg
    assert '@    s3 legend "E(TOT)"' in xvg
    assert "@    s1 hidden true" in xvg
    assert "@    s3 hidden true" not in xvg
    assert "@target G1.S0" in xvg and "@target G2.S0" in xvg

    original = read_energy_table(str(energy_dir / "pq.en"))
    for output in outputs[2:]:
        restored = read_analysis_table(output)
        assert restored.schema.fields == PQ_FIELDS
        np.testing.assert_allclose(restored.data, original.data)
    assert read_analysis_table(outputs[2]).schema == original.schema



def test_convert_energy_xvg_projection_override(energy_dir):
    output = energy_dir / "pq.xvg"

    convert_analysis_output(
        str(energy_dir / "pq.en"),
        [str(output)],
        y_fields=["TEMPERATURE"],
    )

    xvg = output.read_text(encoding="utf-8")
    assert '@    yaxis label "TEMPERATURE [K]"' in xvg
    assert "@    s2 hidden true" in xvg
    assert "@    s1 hidden true" not in xvg
    assert "# XVG_GRAPH TEMPERATURE\n" in xvg
    assert "@g1 on" not in xvg



def test_convert_energy_flag_reads_other_suffix(energy_dir):
    shutil.copy(energy_dir / "pq.en", energy_dir / "run.txt")
    output = energy_dir / "run.csv"

    convert_analysis_output(
        str(energy_dir / "run.txt"),
        [str(output)],
        energy=EnergyFileOptions(
            as_energy_file=True,
            info_file=str(energy_dir / "pq.info"),
        ),
    )

    header = output.read_text(encoding="utf-8").splitlines()[0]
    assert header == PQ_CSV_HEADER



def test_energy_options_need_energy_input(energy_dir):
    with pytest.raises(AnalysisOutputError, match="energy-file input"):
        convert_analysis_output(
            str(energy_dir / "pq.info"),
            [str(energy_dir / "out.csv")],
            energy=EnergyFileOptions(timestep_fs=1.0),
        )



def test_energy_table_rejects_names_equal_after_whitespace(energy_dir):
    (energy_dir / "dup.en").write_text("1 2\n3 4\n", encoding="utf-8")
    (energy_dir / "dup.info").write_text(
        "header\nheader\nheader\n"
        "| E QM 1.0 E_QM 2.0 |\n"
        "-----------------------------------------\n\n",
        encoding="utf-8",
    )

    with pytest.raises(AnalysisOutputError, match="duplicate column names"):
        read_energy_table(
            str(energy_dir / "dup.en"),
            EnergyFileOptions(engine_format="qmcfc"),
        )



def _energy_schema(columns):
    return AnalysisSchema(
        title="run",
        columns=tuple(
            AnalysisColumn(field, field, unit) for field, unit in columns
        ),
        units_in_header=True,
    )



def test_default_plot_skips_zero_energies_and_caps_graphs():
    schema = _energy_schema(
        [("SIMULATION-TIME", "ps"), ("TEMPERATURE", "K")] +
        [(f"E({name})", "kcal/mol") for name in "ABCDEF"]
    )
    data = np.array(
        [
            [1.0, 300.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            [2.0, 301.0, 0.0, 1.5, 2.5, 3.5, 4.5, 5.5]
        ]
    )

    plot = default_energy_plot(schema, data)

    assert plot.y_fields == ("E(B)", "E(C)", "E(D)", "E(E)")
    assert len(plot.graphs) == MAX_DEFAULT_GRAPHS
    assert plot.x_label == "Simulation time [ps]"



def test_default_plot_without_energies_groups_by_unit():
    schema = _energy_schema(
        [
            ("STEP", "1"), ("TEMPERATURE", "K"), ("PRESSURE", "bar"),
            ("T2", "K")
        ]
    )
    data = np.array([[1.0, 300.0, 1.0, 2.0]])

    plot = default_energy_plot(schema, data)

    assert plot.x_field == "STEP"
    assert plot.x_label == "Simulation time"
    assert [graph.y_fields for graph in plot.graphs] == [
        ("TEMPERATURE", "T2"), ("PRESSURE", )
    ]
    assert [graph.y_label
            for graph in plot.graphs] == ["Temperature [K]", "PRESSURE [bar]"]



def test_plot_preset_labels_for_shared_and_unknown_units():
    schema = _energy_schema(
        [
            ("t", "ps"), ("E(TOT)", "kcal/mol"), ("E(KIN)", "kcal/mol"),
            ("VOLUME", "A^3"), ("X", "m"), ("Y", "m"), ("E(QM)", "1"),
            ("QM_MOLECULES", "1")
        ]
    )

    plot = energy_plot_preset(
        schema,
        "VOLUME", ("E(TOT)", "E(KIN)", "X", "Y", "E(QM)", "QM_MOLECULES")
    )

    assert plot.x_label == r"VOLUME [\cE\C\S3\N]"
    assert [graph.y_label for graph in plot.graphs
            ] == ["Energy [kcal/mol]", "Value [m]", "E(QM)", "QM_MOLECULES"]
    assert plot.y_label == "Energy [kcal/mol]"
    assert energy_plot_preset(schema, None, ("X", )).x_label == "Row"



def test_mixed_unit_selection_writes_stacked_graphs(energy_dir):
    output = energy_dir / "pq.xvg"

    convert_analysis_output(
        str(energy_dir / "pq.en"),
        [str(output)],
        y_fields=["E(TOT)", "TEMPERATURE", "E(KIN)"],
    )

    xvg = output.read_text(encoding="utf-8")
    assert "# XVG_GRAPH E(TOT) E(KIN)\n# XVG_GRAPH TEMPERATURE\n" in xvg
    assert "@g1 on" in xvg and "@g2 on" not in xvg
    assert '@    yaxis label "Energy [kcal/mol]"' in xvg
    assert '@    yaxis label "TEMPERATURE [K]"' in xvg
    assert xvg.count("@target G1.") == 1
    assert xvg.count("@target G0.") == len(PQ_FIELDS) - 1
    assert xvg.count("@autoscale") == 2

    restored = read_analysis_table(str(output))
    assert restored.schema.plot.y_fields == ("E(TOT)", "TEMPERATURE", "E(KIN)")
    assert [graph.y_fields for graph in restored.schema.plot.graphs] == [
        ("E(TOT)", "E(KIN)"), ("TEMPERATURE", )
    ]
    original = read_energy_table(str(energy_dir / "pq.en"))
    np.testing.assert_allclose(restored.data, original.data)
    assert restored.schema.columns == original.schema.columns



def test_energy_xvg_converts_back_to_energy_xvg(energy_dir):
    first, second = energy_dir / "a.xvg", energy_dir / "b.xvg"
    convert_analysis_output(str(energy_dir / "pq.en"), [str(first)])

    convert_analysis_output(str(first), [str(second)])

    assert first.read_text(encoding="utf-8"
                           ) == second.read_text(encoding="utf-8")



def test_energy_csv_header_round_trip(energy_dir):
    output = energy_dir / "pq.csv"
    convert_analysis_output(str(energy_dir / "pq.en"), [str(output)])

    restored = read_analysis_table(str(output))

    original = read_energy_table(str(energy_dir / "pq.en"))
    assert restored.schema.columns == tuple(
        AnalysisColumn(column.field, column.field, column.unit)
        for column in original.schema.columns
    )
    assert restored.schema.units_in_header
    rewritten = energy_dir / "again.csv"
    convert_analysis_output(str(output), [str(rewritten)])
    assert rewritten.read_text(encoding="utf-8"
                               ).splitlines()[0] == (PQ_CSV_HEADER)



def test_csv_header_columns_are_pandas_style_names(energy_dir):
    output = energy_dir / "pq.csv"
    convert_analysis_output(str(energy_dir / "pq.en"), [str(output)])

    with output.open(encoding="utf-8", newline="") as file:
        header = next(csv.reader(file))

    assert header[3] == "E(TOT) [kcal/mol]"
    assert header[5] == "N(QM-ATOMS)"
    assert all(name == name.strip() for name in header)



@pytest.mark.parametrize(
    ("unit", "expected"),
    (
        ("A^3", r"\cE\C\S3\N"),
        ("g/cm^3", r"g/cm\S3\N"),
        ("amuA/fs", "amuA/fs"),
        ("kcal/mol", "kcal/mol"),
    ),
)
def test_grace_unit(unit, expected):
    assert grace_unit(unit) == expected



@pytest.mark.parametrize(
    ("values", "expected"),
    (
        ([-186000.0, -185990.0], 0),
        ([4.2, 4.8], 1),
        ([0.8345, 0.8355], 4),
        ([1.0, 1.0], None),
        ([1e9, 2e9], None),
        ([1e-12, 2e-12], None),
        ([np.nan], None),
    ),
)
def test_tick_precision_follows_data_range(values, expected):
    assert _tick_precision(np.array(values)[:, None]) == expected



def test_default_xvg_uses_decimal_ticks_for_large_energies(energy_dir):
    output = energy_dir / "pq.xvg"
    convert_analysis_output(str(energy_dir / "pq.en"), [str(output)])

    xvg = output.read_text(encoding="utf-8")

    assert "@    yaxis ticklabel format decimal" in xvg
    assert "@    xaxis ticklabel off" in xvg
    assert '@    title "PQ energy file pq"' in xvg



def test_write_xvg_rejects_missing_graph_labels():
    document = XVGDocument(
        title="energy",
        fields=("TIME", "E"),
        symbols=("t", "E"),
        units=("ps", "kcal/mol"),
        x_field="TIME",
        y_fields=("E", ),
        data=np.array([[0.0, 1.0], [1.0, 2.0]]),
        graphs=(("E", ), ),
    )
    style = XVGPlotStyle(
        x_label="t",
        y_label="E",
        legends=("t", "E"),
        x_values=np.array([0.0, 1.0]),
    )

    with pytest.raises(XVGFormatError, match="graph labels"):
        write_xvg(io.StringIO(), document, style)
