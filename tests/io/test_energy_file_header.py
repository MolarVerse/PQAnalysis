"""
Tests for the column header of PQ energy files.
"""

import logging
import shutil
from pathlib import Path

import numpy as np
import pytest

from PQAnalysis import __package_name__
from PQAnalysis.analysis.output import convert_analysis_output
from PQAnalysis.io import EnergyFileReader
from PQAnalysis.io.energy_file_header import read_energy_header
from PQAnalysis.io.exceptions import EnergyFileReaderError

from ..conftest import caplog_for_logger
from . import pytestmark  # pylint: disable=unused-import

DATA_DIR = Path(__file__).parents[1] / "data" / "readEnergyFile"
HEADER = "# PQ energy file v1.1.0\n"
LOGGER = f"{__package_name__}.EnergyFileReader"
NAMES = (
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
    for name in (
        "md-01.en", "md-01.info", "header_only.en", "header_timestep.en"
    ):
        shutil.copy(DATA_DIR / name, tmp_path / name)
    return tmp_path



def logged(caplog):
    """Return the captured log text without line wrapping."""
    return " ".join(caplog.text.split())



def read(path, **kwargs):
    """Read one energy file."""
    return EnergyFileReader(str(path), **kwargs).read()



def test_header_labels(energy_dir):
    info, units = read_energy_header(str(energy_dir / "header_only.en"))

    assert tuple(info) == NAMES
    assert list(info.values()) == list(range(12))
    assert units["SIMULATION-TIME"] == "step"
    assert units["N(QM-ATOMS)"] == "-"
    assert units["DENSITY"] == "g/cm^3"



@pytest.mark.parametrize("timestep_fs", (None, 0.5))
def test_header_without_info_matches_info_path(energy_dir, timestep_fs):
    energy = read(energy_dir / "header_only.en", timestep_fs=timestep_fs)
    legacy = read(energy_dir / "md-01.en", timestep_fs=timestep_fs)

    assert energy.info == legacy.info
    assert energy.units == legacy.units
    np.testing.assert_array_equal(energy.data, legacy.data)



def test_header_without_timestep_keeps_steps(energy_dir):
    energy = read(energy_dir / "header_only.en")

    assert energy.simulation_time_unit == "step"
    np.testing.assert_array_equal(energy.simulation_time, [1, 2])
    assert energy.temperature_unit == "K"



def test_timestep_line_inside_header_gives_ps(energy_dir):
    energy = read(energy_dir / "header_timestep.en")

    assert tuple(energy.info) == NAMES
    assert energy.simulation_time_unit == "ps"
    np.testing.assert_allclose(energy.simulation_time, [0.0005, 0.001])



def test_time_column_in_ps_is_not_scaled(tmp_path):
    filename = tmp_path / "run.en"
    filename.write_text(
        HEADER + "#SIMULATION-TIME[ps] TEMPERATURE[K]\n0.5 300\n1.0 301\n",
        encoding="utf-8",
    )

    energy = read(filename, timestep_fs=2.0)

    assert energy.simulation_time_unit == "ps"
    np.testing.assert_array_equal(energy.simulation_time, [0.5, 1.0])



def test_header_agreeing_with_info_does_not_warn(energy_dir, caplog):
    shutil.copy(energy_dir / "md-01.info", energy_dir / "header_only.info")

    with caplog_for_logger(caplog, LOGGER, logging.WARNING):
        energy = read(energy_dir / "header_only.en")

    assert "differ from the info file" not in logged(caplog)
    assert tuple(energy.info) == NAMES



def test_header_wins_over_disagreeing_info(energy_dir, caplog):
    (energy_dir / "run.en").write_text(
        HEADER + "#SIMULATION-TIME[step] TEMPERATURE[K] E(TOT)[kJ/mol]\n"
        "1 300 -5\n",
        encoding="utf-8",
    )
    (energy_dir / "run.info").write_text(
        "header\nheader\nheader\n"
        "| SIMULATION-TIME 0.1 ps TEMPERATURE 300 K |\n"
        "| E(TOT) -5 kcal/mol |\n"
        "-----------------------------------------\n",
        encoding="utf-8",
    )

    with caplog_for_logger(caplog, LOGGER, logging.WARNING):
        energy = read(energy_dir / "run.en")

    assert energy.units["E(TOT)"] == "kJ/mol"
    assert "differ from the info file" in logged(caplog)



def test_header_differs_from_info_by_name(energy_dir, caplog):
    shutil.copy(energy_dir / "md-01.info", energy_dir / "run.info")
    (energy_dir / "run.en").write_text(
        HEADER + "#" + " ".join(f"C{index}[-]" for index in range(12)) + "\n" +
        " ".join("1" for _ in range(12)) + "\n",
        encoding="utf-8",
    )

    with caplog_for_logger(caplog, LOGGER, logging.WARNING):
        energy = read(energy_dir / "run.en")

    assert "C0" in energy.info
    assert "differ from the info file" in logged(caplog)



def test_file_without_header(energy_dir):
    assert read_energy_header(str(energy_dir / "md-01.en")) is None
    assert read(energy_dir / "md-01.en").simulation_time_unit == "step"



def test_multiple_files_share_header(energy_dir):
    energy = EnergyFileReader(
        [
            str(energy_dir / "header_only.en"),
            str(energy_dir / "header_only.en"),
        ]
    ).read()

    np.testing.assert_array_equal(energy.simulation_time, [1, 2, 1, 2])



def test_multiple_files_with_different_headers(energy_dir):
    (energy_dir / "other.en").write_text(
        HEADER + "#" + " ".join(f"C{index}[-]" for index in range(12)) + "\n" +
        " ".join("1" for _ in range(12)) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(EnergyFileReaderError, match="different column"):
        EnergyFileReader(
            [
                str(energy_dir / "header_only.en"),
                str(energy_dir / "other.en"),
            ]
        ).read()



@pytest.mark.parametrize(
    ("column_lines", "message"),
    (
        ("", "needs a column line"),
        ("#\n", "needs a column line"),
        ("# timestep = 0.5 fs\n", "Invalid column 'timestep'"),
        ("#STEP[-] TEMPERATURE\n", "Invalid column 'TEMPERATURE'"),
        ("#STEP[-] T[K]x\n", r"Invalid column 'T\[K\]x'"),
        ("#STEP[-] STEP[-]\n", "Duplicate column 'STEP'"),
        (
            "#STEP[-] SIMULATION-TIME[step]\n",
            "Duplicate column 'SIMULATION-TIME'",
        ),
    ),
)
def test_malformed_header(tmp_path, column_lines, message):
    filename = tmp_path / "bad.en"
    filename.write_text(HEADER + column_lines + "1 2\n", encoding="utf-8")

    with pytest.raises(EnergyFileReaderError, match=message):
        read_energy_header(str(filename))



def test_header_column_count_must_match_data(tmp_path):
    filename = tmp_path / "short.en"
    filename.write_text(
        HEADER + "#SIMULATION-TIME[step] T[K]\n1 2 3\n", encoding="utf-8"
    )

    with pytest.raises(EnergyFileReaderError, match="names 2 columns"):
        read(filename)



@pytest.mark.parametrize(
    "lines",
    (
        "# PQ energy file\n#A[-]\n",
        "#PQ energy file 2\n\n# note\n#A[-]\n\n",
    ),
)
def test_header_line_variants(tmp_path, lines):
    filename = tmp_path / "run.en"
    filename.write_text(lines + "1\n", encoding="utf-8")

    assert read_energy_header(str(filename)) == ({"A": 0}, {"A": "-"})



def test_step_with_other_unit_keeps_its_name(tmp_path):
    filename = tmp_path / "run.en"
    filename.write_text(HEADER + "#STEP[1] T[K]\n1 2\n", encoding="utf-8")

    assert read_energy_header(str(filename)) == (
        {"STEP": 0, "T": 1}, {"STEP": "1", "T": "K"}
    )



def test_other_comment_is_not_a_header(tmp_path):
    filename = tmp_path / "run.en"
    filename.write_text("# PQ energy files\n#A[-]\n1\n", encoding="utf-8")

    assert read_energy_header(str(filename)) is None



def test_convert_uses_header_names(energy_dir):
    output_csv = energy_dir / "run.csv"
    output_xvg = energy_dir / "run.xvg"

    convert_analysis_output(
        str(energy_dir / "header_timestep.en"),
        [str(output_csv), str(output_xvg)],
    )

    csv_lines = output_csv.read_text(encoding="utf-8").splitlines()
    assert csv_lines[0].startswith(
        "SIMULATION-TIME [ps],TEMPERATURE [K],PRESSURE [bar],"
    )
    assert csv_lines[0].count(",") == len(NAMES) - 1
    assert csv_lines[1].startswith("0.0005,298.066019722598,")
    xvg = output_xvg.read_text(encoding="utf-8")
    assert '@    xaxis label "Simulation time [ps]"' in xvg
    assert "@    s0 hidden true" in xvg
    assert '@    s1 legend "TEMPERATURE [K]"' in xvg



def step_header_copy(energy_dir, name):
    """Copy a header fixture with the first column named STEP[-]."""
    text = (energy_dir / name).read_text(encoding="utf-8")
    filename = energy_dir / f"step_{name}"
    filename.write_text(
        text.replace("#SIMULATION-TIME[step] ", "#STEP[-] ", 1),
        encoding="utf-8",
    )
    return filename



def test_step_header_labels(energy_dir):
    filename = step_header_copy(energy_dir, "header_only.en")

    assert "#STEP[-] " in filename.read_text(encoding="utf-8")
    assert read_energy_header(str(filename)) == read_energy_header(
        str(energy_dir / "header_only.en")
    )



@pytest.mark.parametrize(
    ("name", "timestep_fs"),
    (
        ("header_only.en", None),
        ("header_only.en", 0.5),
        ("header_timestep.en", None),
    ),
)
def test_step_header_reads_like_old_header(energy_dir, name, timestep_fs):
    energy = read(
        step_header_copy(energy_dir, name), timestep_fs=timestep_fs
    )
    old = read(energy_dir / name, timestep_fs=timestep_fs)

    assert energy.info == old.info
    assert energy.units == old.units
    assert energy.simulation_time_unit == old.simulation_time_unit
    np.testing.assert_array_equal(energy.data, old.data)



def test_step_header_agreeing_with_info_does_not_warn(energy_dir, caplog):
    filename = step_header_copy(energy_dir, "header_only.en")
    shutil.copy(energy_dir / "md-01.info", filename.with_suffix(".info"))

    with caplog_for_logger(caplog, LOGGER, logging.WARNING):
        energy = read(filename)

    assert "differ from the info file" not in logged(caplog)
    assert tuple(energy.info) == NAMES



def test_convert_step_header_like_old_header(energy_dir):
    filename = step_header_copy(energy_dir, "header_timestep.en")
    new_csv = energy_dir / "new.csv"
    old_csv = energy_dir / "old.csv"

    convert_analysis_output(str(filename), [str(new_csv)])
    convert_analysis_output(
        str(energy_dir / "header_timestep.en"), [str(old_csv)]
    )

    assert new_csv.read_text(encoding="utf-8") == old_csv.read_text(
        encoding="utf-8"
    )
