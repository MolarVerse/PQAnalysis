import pytest
import numpy as np

from . import pytestmark

from collections import defaultdict

from PQAnalysis.io import EnergyFileReader, InfoFileReader
from PQAnalysis.io.exceptions import EnergyFileReaderError
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.traj.exceptions import MDEngineFormatError
from PQAnalysis.exceptions import PQFileNotFoundError



class TestEnergyReader:

    @pytest.mark.parametrize("example_dir", ["readEnergyFile"], indirect=False)
    def test__init__(self, test_with_data_dir):
        with pytest.raises(PQFileNotFoundError) as exception:
            EnergyFileReader("tmp")
        assert str(exception.value) == "File tmp not found."

        with pytest.raises(PQFileNotFoundError) as exception:
            EnergyFileReader(["md-01.en", "tmp"])
        assert str(exception.value) == (
            "At least one of the given files does not exist. File tmp not found."
        )

        reader = EnergyFileReader("md-01.en")
        assert reader.filename == "md-01.en"
        assert reader.multiple_files == False
        assert reader.info_filename == "md-01.info"
        assert reader.with_info_file == True
        assert reader.format == MDEngineFormat.PQ

        reader = EnergyFileReader(["md-01.en", "md-01_noinfo.en"])
        assert reader.filenames == ["md-01.en", "md-01_noinfo.en"]
        assert reader.multiple_files == True
        assert reader.info_filename == "md-01.info"
        assert reader.with_info_file == True
        assert reader.format == MDEngineFormat.PQ

        reader = EnergyFileReader("md-01.en", use_info_file=False)
        assert reader.filename == "md-01.en"
        assert reader.multiple_files == False
        assert reader.info_filename == None
        assert reader.with_info_file == False
        assert reader.format == MDEngineFormat.PQ

        reader = EnergyFileReader("md-01_noinfo.en")
        assert reader.filename == "md-01_noinfo.en"
        assert reader.multiple_files == False
        assert reader.info_filename == None
        assert reader.with_info_file == False
        assert reader.format == MDEngineFormat.PQ

        with pytest.raises(PQFileNotFoundError) as exception:
            EnergyFileReader(
                "md-01_noinfo.en",
                info_filename="md-01_noinfo.info",
                use_info_file=True
            )
        assert str(exception.value) == "Info File md-01_noinfo.info not found."

        reader = EnergyFileReader(
            "md-01_noinfo.en",
            info_filename="md-01.info"
        )
        assert reader.filename == "md-01_noinfo.en"
        assert reader.info_filename == "md-01.info"
        assert reader.with_info_file == True
        assert reader.format == MDEngineFormat.PQ

        reader = EnergyFileReader("md-01.en", engine_format="qmcfc")
        assert reader.filename == "md-01.en"
        assert reader.info_filename == "md-01.info"
        assert reader.with_info_file == True
        assert reader.format == MDEngineFormat.QMCFC

        with pytest.raises(MDEngineFormatError) as exception:
            EnergyFileReader("md-01.en", engine_format="tmp")
        assert str(exception.value) == (
            "\n"
            "'tmp' is not a valid MDEngineFormat.\n"
            f"Possible values are: {MDEngineFormat.member_repr()} "
            "or their case insensitive string representation: "
            f"{MDEngineFormat.value_repr()}"
        )

    @pytest.mark.parametrize("example_dir", ["readEnergyFile"], indirect=False)
    def test__info_file_found__(self, test_with_data_dir):
        reader = EnergyFileReader("md-01.en")
        assert reader.__info_file_found__() == True

        reader = EnergyFileReader("md-01_noinfo.en")
        assert reader.__info_file_found__() == False

        reader = EnergyFileReader(
            "md-01_noinfo.en",
            info_filename="md-01.info"
        )
        assert reader.__info_file_found__() == True

        with pytest.raises(PQFileNotFoundError) as exception:
            EnergyFileReader(
                "md-01_noinfo.en",
                info_filename="md-01_noinfo.info",
                use_info_file=True
            )
        assert str(exception.value) == "Info File md-01_noinfo.info not found."

    @pytest.mark.parametrize("example_dir", ["readEnergyFile"], indirect=False)
    def test_read(self, test_with_data_dir):
        data_ref = np.array(
            [
            [1.00000000e+00,
            2.00000000e+00],
            [2.98066020e+02,
            2.98442926e+02],
            [2.32566666e+04,
            2.39839997e+04],
            [-1.85898221e+05,
            -1.85719252e+05],
            [-1.85903822e+05,
            -1.86024132e+05],
            [0.00000000e+00,
            0.00000000e+00],
            [5.60049937e+00,
            3.04879280e+02],
            [0.00000000e+00,
            0.00000000e+00],
            [3.72913867e+03,
            3.72146189e+03],
            [8.34545285e-01,
            8.34545285e-01],
            [1.12687000e-16,
            5.63381000e-17],
            [1.07571000e+00,
            9.97330000e-01]
            ]
        )

        infoFileReader = InfoFileReader("md-01.info")
        info, units = infoFileReader.read()

        reader = EnergyFileReader("md-01.en")
        energy = reader.read()

        assert energy.data.shape == (12, 2)
        assert np.allclose(energy.data, data_ref)
        assert energy.info == info
        assert energy.units == units
        assert energy.info_given == True
        assert energy.units_given == True

        reader = EnergyFileReader("qmmm.en")
        energy = reader.read()
        assert energy.data.shape == (13, 1)
        assert len(energy.info) == 13
        assert energy.info["N(SM-MOL)"] == 10
        assert energy.units["N(SM-MOL)"] == "-"
        assert energy.data[energy.info["N(SM-MOL)"], 0] == 0
        assert energy.info["MOMENTUM"] == 11
        assert energy.info["LOOPTIME"] == 12

        reader = EnergyFileReader("md-01_noinfo.en")
        energy = reader.read()
        assert np.allclose(energy.data, data_ref)
        assert energy.info == defaultdict(lambda: None)
        assert energy.units == defaultdict(lambda: None)
        assert energy.info_given == False
        assert energy.units_given == False

        reader = EnergyFileReader(["md-01.en", "md-01.en"])
        energy = reader.read()
        assert energy.data.shape == (12, 4)
        assert np.allclose(energy.data, np.hstack([data_ref, data_ref]))
        assert energy.info == info
        assert energy.units == units
        assert energy.info_given == True
        assert energy.units_given == True

        reader = EnergyFileReader(
            ["md-01_noinfo.en", "md-01_noinfo.en"],
            use_info_file=False,
        )
        energy = reader.read()
        assert energy.data.shape == (12, 4)
        assert np.allclose(energy.data, np.hstack([data_ref, data_ref]))
        assert energy.info == defaultdict(lambda: None)
        assert energy.units == defaultdict(lambda: None)
        assert energy.info_given == False
        assert energy.units_given == False

    @pytest.mark.usefixtures("tmpdir")
    def test_read_skips_blank_lines(self):
        with open("md-blank.en", "w", encoding="utf-8") as file:
            print("1 2.0 3.0", file=file)
            print("", file=file)
            print("2 4.0 5.0", file=file)
            print("", file=file)

        reader = EnergyFileReader("md-blank.en", use_info_file=False)
        energy = reader.read()

        assert energy.data.shape == (3, 2)
        assert np.allclose(
            energy.data,
            np.array([[1.0, 2.0], [2.0, 4.0], [3.0, 5.0]])
        )

    @pytest.mark.usefixtures("tmpdir")
    def test_read_invalid_number_of_columns(self):
        with open("md-ragged.en", "w", encoding="utf-8") as file:
            print("1 2.0 3.0", file=file)
            print("2 4.0 5.0 6.0", file=file)

        reader = EnergyFileReader("md-ragged.en", use_info_file=False)

        with pytest.raises(EnergyFileReaderError) as exception:
            reader.read()
        assert str(exception.value) == (
            "Invalid number of columns in energy file md-ragged.en line 2. "
            "Expected 3 columns."
        )

    @pytest.mark.usefixtures("tmpdir")
    def test_read_inconsistent_columns_across_files(self):
        with open("md-a.en", "w", encoding="utf-8") as file:
            print("1 2.0 3.0", file=file)

        with open("md-b.en", "w", encoding="utf-8") as file:
            print("1 2.0 3.0 4.0", file=file)

        reader = EnergyFileReader(["md-a.en", "md-b.en"], use_info_file=False)

        with pytest.raises(EnergyFileReaderError) as exception:
            reader.read()
        assert str(exception.value) == (
            "Energy file md-b.en has 4 columns, but the previously read "
            "energy files have 3 columns."
        )
