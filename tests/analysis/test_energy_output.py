"""
Tests for converting PQ and QMCFC energy files to analysis tables.
"""

import shutil
from pathlib import Path

import numpy as np
import pytest

from PQAnalysis.analysis.output import (
    AnalysisOutputError,
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
    assert plot.y_fields == PQ_FIELDS[1:]
    assert plot.x_label == "SIMULATION-TIME [step]"
    assert plot.legends[0] == "TEMPERATURE [K]"



def test_pq_energy_timestep_gives_time_in_ps(energy_dir):
    table = read_energy_table(
        str(energy_dir / "pq.en"), EnergyFileOptions(timestep_fs=2.0)
    )

    assert table.schema.columns[0].unit == "ps"
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
    assert table.schema.plot.x_field == "SIMULATION_TIME"
    assert table.schema.plot.legends[0] == "QM_MOLECULES"



def test_energy_table_without_info(energy_dir):
    (energy_dir / "pq.info").unlink()

    table = read_energy_table(str(energy_dir / "pq.en"))

    assert table.schema.fields[:2] == ("column_1", "column_2")
    assert table.schema.plot.x_field == "column_1"
    assert len(table.schema.plot.y_fields) == 11



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
    assert csv_lines[0] == ",".join(PQ_FIELDS)
    assert csv_lines[1].startswith("1.0,298.066019722598,")
    tsv_header = (energy_dir / "pq.tsv").read_text(encoding="utf-8")
    assert tsv_header.startswith("\t".join(PQ_FIELDS) + "\n")

    xvg = (energy_dir / "pq.xvg").read_text(encoding="utf-8")
    assert '@    xaxis label "SIMULATION-TIME [step]"' in xvg
    assert '@    s1 legend "TEMPERATURE [K]"' in xvg
    assert '@    s3 legend "E(TOT) [kcal/mol]"' in xvg
    assert "@    s0 hidden true" in xvg
    assert "@    s1 hidden true" not in xvg

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
    assert header == ",".join(PQ_FIELDS)



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
