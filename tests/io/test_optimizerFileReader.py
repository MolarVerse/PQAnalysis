import numpy as np
import pytest

from . import pytestmark

from PQAnalysis.exceptions import PQFileNotFoundError
from PQAnalysis.io import OptimizerFileReader, read_optimizer_file
from PQAnalysis.io.exceptions import OptimizerReaderError



class TestOptimizerFileReader:

    @pytest.mark.parametrize(
        "example_dir", ["readOptimizerFile"], indirect=False
    )
    def test__init__(self, test_with_data_dir):
        with pytest.raises(PQFileNotFoundError) as exception:
            OptimizerFileReader("missing.opt")
        assert str(exception.value) == "File missing.opt not found."

        reader = OptimizerFileReader("optimization.opt")

        assert reader.filename == "optimization.opt"
        assert reader.multiple_files is False

    @pytest.mark.parametrize(
        "example_dir", ["readOptimizerFile"], indirect=False
    )
    def test_read_real_pq_output(self, test_with_data_dir):
        energy = OptimizerFileReader("optimization.opt").read()

        expected_parameters = (
            "SIMULATION-TIME",
            "ABS-ENERGY-CHANGE",
            "REL-ENERGY-CHANGE",
            "MAX-FORCE",
            "RMS-FORCE",
            "REL-ENERGY-CONV",
            "ABS-ENERGY-CONV",
            "MAX-FORCE-CONV",
            "RMS-FORCE-CONV",
            "REL-ENERGY-LIMIT",
            "ABS-ENERGY-LIMIT",
            "MAX-FORCE-LIMIT",
            "RMS-FORCE-LIMIT",
        )

        assert energy.data.shape == (13, 11)
        assert energy.info == {
            parameter: index
            for index, parameter in enumerate(expected_parameters)
        }
        assert energy.units == {
            "SIMULATION-TIME": "step",
            "ABS-ENERGY-CHANGE": "kcal/mol",
            "REL-ENERGY-CHANGE": "-",
            "MAX-FORCE": "kcal/mol/A",
            "RMS-FORCE": "kcal/mol/A",
            "REL-ENERGY-CONV": "state",
            "ABS-ENERGY-CONV": "state",
            "MAX-FORCE-CONV": "state",
            "RMS-FORCE-CONV": "state",
            "REL-ENERGY-LIMIT": "-",
            "ABS-ENERGY-LIMIT": "kcal/mol",
            "MAX-FORCE-LIMIT": "kcal/mol/A",
            "RMS-FORCE-LIMIT": "kcal/mol/A",
        }
        assert np.array_equal(energy.simulation_time, np.arange(1, 12))
        assert energy.simulation_time_unit == "step"
        assert energy.data[1, -1] == pytest.approx(388.243509)
        assert np.array_equal(energy.data[5:9, -1], np.full(4, -1.0))

    @pytest.mark.usefixtures("tmpdir")
    def test_read_optimizer_file_with_one_data_row(self):
        with open("optimization.opt", "w", encoding="utf-8") as file:
            print("# PQ optimizer output", file=file)
            print("", file=file)
            print("1 2 3 4 5 6 7 8 9 10 11 12 13", file=file)

        energy = read_optimizer_file("optimization.opt")

        assert energy.data.shape == (13, 1)
        assert np.array_equal(energy.data[:, 0], np.arange(1, 14))

    @pytest.mark.usefixtures("tmpdir")
    def test_invalid_number_of_columns(self):
        with open("optimization.opt", "w", encoding="utf-8") as file:
            print("# PQ optimizer output", file=file)
            print("", file=file)
            print("1 2 3", file=file)

        with pytest.raises(OptimizerReaderError) as exception:
            OptimizerFileReader("optimization.opt").read()
        assert str(exception.value) == (
            "Invalid number of columns in optimizer file line 3. "
            "Expected 13 columns."
        )

    @pytest.mark.usefixtures("tmpdir")
    def test_invalid_numeric_value(self):
        with open("optimization.opt", "w", encoding="utf-8") as file:
            print("1 2 3 4 5 1 1 1 1 10 11 bad 13", file=file)

        with pytest.raises(OptimizerReaderError) as exception:
            OptimizerFileReader("optimization.opt").read()
        assert str(exception.value) == (
            "Invalid numeric value in optimizer file line 1: "
            "1 2 3 4 5 1 1 1 1 10 11 bad 13"
        )

    @pytest.mark.usefixtures("tmpdir")
    def test_empty_file(self):
        with open("optimization.opt", "w", encoding="utf-8") as file:
            print("# PQ optimizer output", file=file)

        with pytest.raises(OptimizerReaderError) as exception:
            OptimizerFileReader("optimization.opt").read()
        assert str(exception.value) == (
            "Optimizer file optimization.opt does not contain optimizer data."
        )
