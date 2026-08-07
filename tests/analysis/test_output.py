"""
Tests for shared analysis-table input, output and conversion.
"""

import numpy as np
import pytest

from PQAnalysis.analysis.output import (
    AnalysisDataWriter,
    AnalysisOutputError,
    AnalysisOutputFormat,
    AnalysisTable,
    AnalysisTableWriter,
    MSD_SCHEMA,
    RDF_SCHEMA,
    convert_analysis_output,
    infer_output_format,
    normal_mode_schema,
    read_analysis_table,
)
from PQAnalysis.analysis.rdf.rdf_output_file_writer import RDFDataWriter
from PQAnalysis.io.exceptions import FileWritingModeError



@pytest.fixture
def rdf_columns():
    """Return a small RDF data set in schema order."""
    return (
        np.array([0.25, 0.75]),
        np.array([0.0, 1.5]),
        np.array([0.0, 2.0]),
        np.array([0.0, 3.25]),
        np.array([-1.0, 0.5]),
    )



@pytest.fixture
def rdf_table(rdf_columns):
    """Return a small validated RDF table."""
    return AnalysisTable.from_columns(RDF_SCHEMA, rdf_columns)



@pytest.mark.parametrize(
    ("filename", "expected"),
    (
        (None, AnalysisOutputFormat.NATIVE),
        ("rdf.dat", AnalysisOutputFormat.NATIVE),
        ("rdf.out", AnalysisOutputFormat.NATIVE),
        ("rdf", AnalysisOutputFormat.NATIVE),
        ("rdf.CSV", AnalysisOutputFormat.CSV),
        ("rdf.tsv", AnalysisOutputFormat.TSV),
        ("rdf.xvg", AnalysisOutputFormat.XVG),
    ),
)
def test_infer_output_format(filename, expected):
    assert infer_output_format(filename) == expected



@pytest.mark.parametrize(
    ("suffix", "header"),
    (
        ("csv", ",".join(RDF_SCHEMA.fields)),
        ("tsv", "\t".join(RDF_SCHEMA.fields)),
    ),
)
def test_delimited_round_trip(tmp_path, rdf_table, suffix, header):
    output_file = tmp_path / f"rdf.{suffix}"

    AnalysisTableWriter(str(output_file)).write(rdf_table)

    assert output_file.read_text(encoding="utf-8").splitlines()[0] == header
    restored = read_analysis_table(str(output_file))
    assert restored.schema == RDF_SCHEMA
    np.testing.assert_allclose(restored.data, rdf_table.data)



@pytest.mark.parametrize("suffix", ("csv", "tsv"))
def test_single_column_delimited_round_trip(tmp_path, suffix):
    output_file = tmp_path / f"mode.{suffix}"
    table = AnalysisTable.from_columns(
        normal_mode_schema(1),
        (np.array([0.25, -0.5]), ),
    )

    AnalysisTableWriter(str(output_file)).write(table)

    restored = read_analysis_table(str(output_file))
    assert restored.schema == table.schema
    np.testing.assert_allclose(restored.data, table.data)



def test_nonsequential_mode_fields_remain_generic(tmp_path):
    input_file = tmp_path / "modes.csv"
    input_file.write_text(
        "mode_2,mode_1\n0.25,-0.5\n",
        encoding="utf-8",
    )

    table = read_analysis_table(str(input_file))

    assert table.schema.fields == ("mode_2", "mode_1")



def test_native_round_trip_for_unknown_suffix(tmp_path, rdf_table):
    output_file = tmp_path / "rdf.out"

    AnalysisTableWriter(str(output_file)).write(rdf_table)

    text = output_file.read_text(encoding="utf-8")
    assert text.startswith("# PQAnalysis: Radial distribution function\n")
    assert "# SYMBOLS rᵢ g(rᵢ) N(rᵢ) g(rᵢ)ΔVᵢ Hᵢ−Eᵢ\n" in text
    assert read_analysis_table(str(output_file)).schema == RDF_SCHEMA



def test_reader_detects_csv_content_with_dat_suffix(tmp_path):
    input_file = tmp_path / "table.dat"
    input_file.write_text(
        "r_i,g_r_i,N_r_i,g_r_i_dV_i,H_i_minus_E_i\n"
        "0.25,0.0,0.0,0.0,-1.0\n",
        encoding="utf-8",
    )

    table = read_analysis_table(str(input_file))

    assert table.schema == RDF_SCHEMA
    np.testing.assert_allclose(table.data[0], (0.25, 0.0, 0.0, 0.0, -1.0))



def test_reader_accepts_legacy_headerless_native_data(tmp_path):
    input_file = tmp_path / "legacy.dat"
    input_file.write_text("0.0 1.0\n1.0 2.0\n", encoding="utf-8")

    table = read_analysis_table(str(input_file))

    assert table.schema.fields == ("column_1", "column_2")
    assert table.schema.title == "legacy"
    np.testing.assert_allclose(table.data, ((0.0, 1.0), (1.0, 2.0)))



@pytest.mark.parametrize(
    ("contents", "message"),
    (
        ("", "is empty"),
        ("# comment only\n", "contains no numeric rows"),
        ("a,a\n1,2\n", "duplicate field names"),
        ("a,\n1,2\n", "blank field name"),
        ("a,b\n1\n", "inconsistent row widths"),
        ("a,b\n1,nope\n", "invalid numeric data"),
        (
            "# PQAnalysis: Broken\n"
            "# FIELDS a b\n"
            "# SYMBOLS a b\n"
            "# UNITS 1\n"
            "1 2\n",
            "inconsistent FIELDS",
        ),
    ),
)
def test_reader_rejects_malformed_tables(tmp_path, contents, message):
    input_file = tmp_path / "broken.dat"
    input_file.write_text(contents, encoding="utf-8")

    with pytest.raises(AnalysisOutputError, match=message):
        read_analysis_table(str(input_file))



def test_rdf_xvg_uses_scientific_plot_preset(tmp_path, rdf_table):
    output_file = tmp_path / "rdf.xvg"

    AnalysisTableWriter(str(output_file)).write(rdf_table)

    text = output_file.read_text(encoding="utf-8")
    assert '@    title "Radial distribution function"' in text
    assert '@    xaxis label "r [\\cE\\C]"' in text
    assert '@    yaxis label "g(r)"' in text
    assert "# PQAnalysis-XVG: 1" in text
    assert "# XVG_X_FIELD r_i" in text
    assert "# XVG_Y_FIELDS g_r_i" in text
    assert '@    s1 legend "g(r)"' in text
    assert "@    s0 hidden true" in text
    assert text.count("@type xy") == 5
    assert text.count("\n&\n") == 5
    assert "@target G0.S0\n@type xy\n0.25 0.25\n0.75 0.75\n&" in text
    assert "@target G0.S1\n@type xy\n0.25 0.0\n0.75 1.5\n&" in text

    restored = read_analysis_table(str(output_file))
    assert restored.schema == rdf_table.schema
    np.testing.assert_allclose(restored.data, rdf_table.data)



def test_multiseries_xvg_preserves_all_columns(tmp_path):
    output_file = tmp_path / "msd.xvg"
    table = AnalysisTable.from_columns(
        MSD_SCHEMA,
        (
            np.array([0, 1]),
            np.array([0.0, 1.0]),
            np.array([0.0, 2.0]),
            np.array([0.0, 3.0]),
        ),
    )

    AnalysisTableWriter(str(output_file)).write(table)

    text = output_file.read_text(encoding="utf-8")
    assert text.count("@type xy") == 4
    assert text.count("\n&\n") == 4
    assert '@    s1 legend "x"' in text
    assert '@    s2 legend "y"' in text
    assert '@    s3 legend "z"' in text
    assert "@    s0 hidden true" in text
    assert "@    s1 hidden true" not in text

    restored = read_analysis_table(str(output_file))
    assert restored.schema == table.schema
    np.testing.assert_allclose(restored.data, table.data)



def test_custom_xvg_projection_uses_grace_safe_scientific_labels(
    tmp_path, rdf_table
):
    output_file = tmp_path / "rdf-volume.xvg"

    AnalysisTableWriter(str(output_file)).write(
        rdf_table,
        y_fields=("g_r_i_dV_i", ),
    )

    text = output_file.read_text(encoding="utf-8")
    label = r"g(r\si\N)DeltaV\si\N [\cE\C\S3\N]"
    assert f'@    yaxis label "{label}"' in text
    assert f'@    s3 legend "{label}"' in text
    assert "@    s3 hidden true" not in text
    assert "# SYMBOLS rᵢ g(rᵢ)" in text
    assert all(
        "ᵢ" not in line for line in text.splitlines() if line.startswith("@")
    )

    restored = read_analysis_table(str(output_file))
    assert restored.schema.plot.y_fields == ("g_r_i_dV_i", )
    np.testing.assert_allclose(restored.data, rdf_table.data)



def test_xvg_round_trip_with_synthetic_row_axis(tmp_path):
    output_file = tmp_path / "modes.xvg"
    table = AnalysisTable(
        schema=normal_mode_schema(2),
        data=np.array(((0.25, -0.5), (0.75, 0.5))),
    )

    AnalysisTableWriter(str(output_file)).write(table)

    text = output_file.read_text(encoding="utf-8")
    assert "# XVG_X_FIELD @row" in text
    assert "@target G0.S0\n@type xy\n1 0.25\n2 0.75\n&" in text
    assert "@target G0.S1\n@type xy\n1 -0.5\n2 0.5\n&" in text
    restored = read_analysis_table(str(output_file))
    assert restored.schema == table.schema
    np.testing.assert_allclose(restored.data, table.data)



def test_xvg_is_detected_without_xvg_extension(tmp_path, rdf_table):
    xvg_file = tmp_path / "rdf.xvg"
    disguised_file = tmp_path / "rdf.out"
    AnalysisTableWriter(str(xvg_file)).write(rdf_table)
    xvg_file.rename(disguised_file)

    restored = read_analysis_table(str(disguised_file))

    assert restored.schema == rdf_table.schema
    np.testing.assert_allclose(restored.data, rdf_table.data)



def test_invalid_xvg_projection_does_not_create_output(tmp_path, rdf_table):
    output_file = tmp_path / "rdf.xvg"
    writer = AnalysisTableWriter(str(output_file))

    with pytest.raises(AnalysisOutputError, match="Unknown analysis field"):
        writer.write(rdf_table, y_fields=("missing", ))

    assert not output_file.exists()



def test_duplicate_xvg_y_fields_do_not_create_output(tmp_path, rdf_table):
    output_file = tmp_path / "rdf.xvg"
    writer = AnalysisTableWriter(str(output_file))

    with pytest.raises(AnalysisOutputError, match="must be distinct"):
        writer.write(rdf_table, y_fields=("g_r_i", "g_r_i"))

    assert not output_file.exists()



@pytest.mark.parametrize(
    ("contents", "message"),
    (
        (
            "# Conventional XVG\n@target G0.S0\n@type xy\n0.0 1.0\n&\n",
            "does not contain reversible PQAnalysis table metadata",
        ),
        (
            "# PQAnalysis-XVG: 1\n"
            "# PQAnalysis: Broken\n"
            "# FIELDS x y\n"
            "# SYMBOLS x y\n"
            "# UNITS 1 1\n"
            "# XVG_X_FIELD x\n"
            "# XVG_Y_FIELDS y\n"
            "@target G0.S0\n"
            "@type xy\n"
            "0.0 0.0\n"
            "&\n",
            "inconsistent FIELDS",
        ),
        (
            "# PQAnalysis-XVG: 1\n"
            "# PQAnalysis: Broken\n"
            "# FIELDS x y\n"
            "# SYMBOLS x y\n"
            "# UNITS 1 1\n"
            "# XVG_X_FIELD x\n"
            "# XVG_Y_FIELDS y\n"
            "@target G0.S0\n"
            "@type xy\n"
            "9.0 0.0\n"
            "&\n"
            "@target G0.S1\n"
            "@type xy\n"
            "9.0 1.0\n"
            "&\n",
            "inconsistent x-axis",
        ),
        (
            "# PQAnalysis-XVG: 1\n"
            "# PQAnalysis: Broken\n"
            "# FIELDS x y\n"
            "# SYMBOLS x y\n"
            "# UNITS 1 1\n"
            "# XVG_X_FIELD x\n"
            "# XVG_Y_FIELDS y\n"
            "@target G0.S0\n"
            "@type xy\n"
            "0.0 0.0\n"
            "&\n"
            "@target G0.S1\n"
            "@type xy\n"
            "1.0 1.0\n"
            "&\n",
            "inconsistent x-axis",
        ),
        (
            "# PQAnalysis-XVG: 1\n"
            "# PQAnalysis: Broken\n"
            "# FIELDS x\n"
            "# SYMBOLS x\n"
            "# UNITS 1\n"
            "# XVG_X_FIELD x\n"
            "# XVG_Y_FIELDS x\n"
            "@target G0.S0\n"
            "@type xy\n"
            "0.0 0.0\n",
            "complete XY data-set blocks",
        ),
        (
            "# PQAnalysis-XVG: 1\n"
            "# PQAnalysis: Broken\n"
            "# FIELDS x\n"
            "# SYMBOLS x\n"
            "# UNITS 1\n"
            "# XVG_X_FIELD x\n"
            "# XVG_Y_FIELDS x\n"
            "@target G0.S1\n"
            "@type xy\n"
            "0.0 0.0\n"
            "&\n",
            "sequential Grace sets",
        ),
        (
            "# PQAnalysis-XVG: 1\n"
            "# PQAnalysis: Broken\n"
            "# FIELDS x\n"
            "# SYMBOLS x\n"
            "# UNITS 1\n"
            "# XVG_X_FIELD x\n"
            "# XVG_Y_FIELDS x\n"
            "@target G0.S0\n"
            "@type xydy\n"
            "0.0 0.0\n"
            "&\n",
            "must use XY data sets",
        ),
    ),
)
def test_reader_rejects_nonreversible_or_malformed_xvg(
    tmp_path, contents, message
):
    input_file = tmp_path / "broken.xvg"
    input_file.write_text(contents, encoding="utf-8")

    with pytest.raises(AnalysisOutputError, match=message):
        read_analysis_table(str(input_file))



def test_convert_writes_multiple_outputs(tmp_path, rdf_table):
    input_file = tmp_path / "rdf.dat"
    csv_file = tmp_path / "rdf.csv"
    tsv_file = tmp_path / "rdf.tsv"
    xvg_file = tmp_path / "rdf.xvg"
    AnalysisTableWriter(str(input_file)).write(rdf_table)

    convert_analysis_output(
        str(input_file),
        (str(csv_file), str(tsv_file), str(xvg_file)),
    )

    assert read_analysis_table(str(csv_file)).schema == RDF_SCHEMA
    assert read_analysis_table(str(tsv_file)).schema == RDF_SCHEMA
    assert '@    yaxis label "g(r)"' in xvg_file.read_text(encoding="utf-8")



@pytest.mark.parametrize("input_suffix", ("out", "csv", "tsv", "xvg"))
@pytest.mark.parametrize("output_suffix", ("dat", "csv", "tsv", "xvg"))
def test_convert_supports_every_format_direction(
    tmp_path,
    rdf_table,
    input_suffix,
    output_suffix,
):
    input_file = tmp_path / f"input.{input_suffix}"
    output_file = tmp_path / f"output.{output_suffix}"
    AnalysisTableWriter(str(input_file)).write(rdf_table)

    convert_analysis_output(str(input_file), (str(output_file), ))

    restored = read_analysis_table(str(output_file))
    assert restored.schema == rdf_table.schema
    np.testing.assert_allclose(restored.data, rdf_table.data)



def test_convert_validates_all_outputs_before_writing(tmp_path, rdf_table):
    input_file = tmp_path / "rdf.dat"
    csv_file = tmp_path / "rdf.csv"
    xvg_file = tmp_path / "rdf.xvg"
    AnalysisTableWriter(str(input_file)).write(rdf_table)

    with pytest.raises(AnalysisOutputError, match="Unknown analysis field"):
        convert_analysis_output(
            str(input_file),
            (str(csv_file), str(xvg_file)),
            y_fields=("missing", ),
        )

    assert not csv_file.exists()
    assert not xvg_file.exists()



def test_convert_rejects_existing_output_before_writing(tmp_path, rdf_table):
    input_file = tmp_path / "rdf.dat"
    new_file = tmp_path / "rdf.tsv"
    existing_file = tmp_path / "rdf.csv"
    AnalysisTableWriter(str(input_file)).write(rdf_table)
    existing_file.write_text("existing data\n", encoding="utf-8")

    with pytest.raises(
        FileWritingModeError,
        match=r"File .*rdf\.csv already exists",
    ):
        convert_analysis_output(
            str(input_file),
            (str(new_file), str(existing_file)),
        )

    assert not new_file.exists()
    assert existing_file.read_text(encoding="utf-8") == "existing data\n"



@pytest.mark.parametrize("duplicate_output", (False, True))
def test_convert_rejects_colliding_paths(
    tmp_path, rdf_table, duplicate_output
):
    input_file = tmp_path / "rdf.dat"
    output_file = tmp_path / "rdf.csv"
    AnalysisTableWriter(str(input_file)).write(rdf_table)
    outputs = (
        (str(output_file), str(output_file)) if duplicate_output else
        (str(input_file), )
    )

    with pytest.raises(AnalysisOutputError, match="must be distinct"):
        convert_analysis_output(str(input_file), outputs)



def test_append_is_rejected_for_structured_formats(tmp_path):
    with pytest.raises(
        AnalysisOutputError, match="Appending is not supported"
    ):
        AnalysisTableWriter(str(tmp_path / "rdf.csv"), mode="a")

    with pytest.raises(
        AnalysisOutputError, match="Appending is not supported"
    ):
        AnalysisDataWriter(str(tmp_path / "rdf.tsv"), mode="a")



def test_analysis_writer_emits_native_and_multiple_exports(
    tmp_path,
    rdf_columns,
):
    native_file = tmp_path / "rdf.dat"
    csv_file = tmp_path / "rdf.csv"
    tsv_file = tmp_path / "rdf.tsv"
    xvg_file = tmp_path / "rdf.xvg"

    RDFDataWriter(
        str(native_file),
        export_files=(str(csv_file), str(tsv_file), str(xvg_file)),
    ).write(rdf_columns)

    assert native_file.read_text(encoding="utf-8").startswith("# PQAnalysis:")
    assert csv_file.read_text(encoding="utf-8").startswith("r_i,g_r_i,")
    assert tsv_file.read_text(encoding="utf-8").startswith("r_i\tg_r_i\t")
    assert '@    title "Radial distribution function"' in xvg_file.read_text(
        encoding="utf-8"
    )



def test_analysis_writer_uses_primary_file_suffix(tmp_path, rdf_columns):
    output_file = tmp_path / "rdf.csv"

    RDFDataWriter(str(output_file)).write(rdf_columns)

    assert output_file.read_text(encoding="utf-8").startswith("r_i,g_r_i,")
