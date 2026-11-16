"""PQ step columns require explicit timestep metadata before claiming ps."""

import numpy as np
import pytest

from PQAnalysis.io import EnergyFileReader
from PQAnalysis.io.exceptions import EnergyFileReaderError



def write_energy(tmp_path, name="run", header=""):
    """Write two sampled steps, not two consecutive integration steps."""
    energy = tmp_path / f"{name}.en"
    energy.write_text(header + "10 298\n100 300\n", encoding="utf-8")
    (tmp_path / f"{name}.info").write_text(
        "header\nheader\nheader\n"
        "| SIMULATION-TIME 0.05 ps TEMPERATURE 300 K |\n"
        "-----------------------------------------\n",
        encoding="utf-8",
    )
    return str(energy)



def test_legacy_pq_numbers_are_steps_not_summary_time(tmp_path):
    energy = EnergyFileReader(write_energy(tmp_path)).read()
    np.testing.assert_array_equal(energy.simulation_time, [10, 100])
    assert energy.simulation_time_unit == "step"
    assert energy.units["TEMPERATURE"] == "K"
    np.testing.assert_array_equal(energy.temperature, [298, 300])



@pytest.mark.parametrize(
    "header",
    [
        "# timestep = 0.5 fs\n",
        "\n# timestep = 0.5 fs\n\n",
        "# timestep = 5e-4 ps\n",
        "  # timestep = 0.5 fs\n# unrelated comment\n",
        "# timestep = 0.5 fs\n# timestep = 0.0005 ps\n",
    ]
)
def test_header_converts_sampled_step_counts_to_ps(tmp_path, header):
    energy = EnergyFileReader(write_energy(tmp_path, header=header)).read()
    np.testing.assert_allclose(energy.simulation_time, [0.005, 0.05])
    assert energy.simulation_time_unit == "ps"
    np.testing.assert_array_equal(energy.temperature, [298, 300])



def test_explicit_timestep_converts_legacy_file(tmp_path):
    energy = EnergyFileReader(write_energy(tmp_path), timestep_fs=0.5).read()
    np.testing.assert_allclose(energy.simulation_time, [0.005, 0.05])
    assert energy.simulation_time_unit == "ps"



@pytest.mark.parametrize("value", [0, -0.5, float("nan"), float("inf"), True])
def test_invalid_explicit_timestep(tmp_path, value):
    with pytest.raises(EnergyFileReaderError, match="finite and positive"):
        EnergyFileReader(write_energy(tmp_path), timestep_fs=value)



@pytest.mark.parametrize(
    "header, message",
    [
        ("# timestep = 0 fs\n", "finite and positive"),
        ("# timestep = -1 fs\n", "finite and positive"),
        ("# timestep = nan fs\n", "finite and positive"),
        ("# timestep = inf fs\n", "finite and positive"),
        ("# timestep = 1e308 ps\n", "finite and positive"),
        ("# timestep = wrong fs\n", "Invalid timestep"),
        ("# timestep = 0.5\n", "Invalid timestep"),
        ("# timestep = 0.5 seconds\n", "Unsupported timestep unit"),
        ("# timestep = 0.5 fs\n# timestep = 1 fs\n", "Conflicting timestep"),
    ]
)
def test_invalid_header_is_not_silently_ignored(tmp_path, header, message):
    with pytest.raises(EnergyFileReaderError, match=message):
        EnergyFileReader(write_energy(tmp_path, header=header)).read()



def test_explicit_and_header_timestep_must_agree(tmp_path):
    filename = write_energy(tmp_path, header="# timestep = 0.5 fs\n")
    with pytest.raises(
        EnergyFileReaderError, match="conflicts with the header"
    ):
        EnergyFileReader(filename, timestep_fs=1).read()
    energy = EnergyFileReader(filename, timestep_fs=0.5).read()
    np.testing.assert_allclose(energy.simulation_time, [0.005, 0.05])



def test_chained_step_counts_keep_their_offsets(tmp_path):
    first = write_energy(tmp_path, "first", "# timestep = 0.5 fs\n")
    second = write_energy(tmp_path, "second", "# timestep = 0.5 fs\n")
    (tmp_path / "second.en"
     ).write_text("# timestep = 0.5 fs\n110 301\n200 302\n", encoding="utf-8")
    energy = EnergyFileReader([first, second]).read()
    np.testing.assert_allclose(
        energy.simulation_time, [0.005, 0.05, 0.055, 0.1]
    )
    np.testing.assert_array_equal(energy.temperature, [298, 300, 301, 302])



def test_mixed_known_unknown_metadata_requires_explicit_timestep(tmp_path):
    first = write_energy(tmp_path, "first", "# timestep = 0.5 fs\n")
    second = write_energy(tmp_path, "second")
    with pytest.raises(EnergyFileReaderError, match="with and without"):
        EnergyFileReader([first, second]).read()
    energy = EnergyFileReader([first, second], timestep_fs=0.5).read()
    np.testing.assert_allclose(energy.simulation_time, [0.005, 0.05] * 2)



def test_raw_numeric_read_without_info_keeps_columns(tmp_path):
    filename = write_energy(tmp_path, header="# timestep = 0.5 fs\n")
    energy = EnergyFileReader(filename, use_info_file=False).read()
    np.testing.assert_array_equal(energy.data, [[10, 100], [298, 300]])
    assert not energy.info_given
    assert not energy.units_given



def test_qmcfc_numbers_are_not_rescaled_as_pq_steps(tmp_path):
    filename = write_energy(tmp_path, header="# timestep = 0.5 fs\n")
    (tmp_path / "run.info").write_text(
        "header\nheader\nheader\n"
        "| SIMULATION TIME 5 TEMPERATURE 300 |\nfooter\nfooter\n",
        encoding="utf-8",
    )
    energy = EnergyFileReader(filename, engine_format="qmcfc").read()
    np.testing.assert_array_equal(energy.simulation_time, [10, 100])
    assert not energy.units_given
    with pytest.raises(EnergyFileReaderError, match="only supported for PQ"):
        EnergyFileReader(filename, engine_format="qmcfc", timestep_fs=0.5)
