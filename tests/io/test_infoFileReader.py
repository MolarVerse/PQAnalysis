import pytest

from PQAnalysis.io.infoFileReader import InfoFileReader


@pytest.mark.parametrize("example_dir", ["readInfoFile"], indirect=False)
def test__init__(test_with_data_dir):
    with pytest.raises(FileNotFoundError) as exception:
        InfoFileReader("tmp")
    assert str(exception.value) == "File tmp not found."

    assert InfoFileReader("md-01.info").filename == "md-01.info"


@pytest.mark.parametrize("example_dir", ["readInfoFile"], indirect=False)
def test_read(test_with_data_dir):
    reader = InfoFileReader("md-01.info")
    info, units = reader.read()

    assert info["SIMULATION-TIME"] == 0
    assert units["SIMULATION-TIME"] == "ps"

    assert info["TEMPERATURE"] == 1
    assert units["TEMPERATURE"] == "K"

    assert info["PRESSURE"] == 2
    assert units["PRESSURE"] == "bar"

    assert info["E(TOT)"] == 3
    assert units["E(TOT)"] == "kcal/mol"

    assert info["E(QM)"] == 4
    assert units["E(QM)"] == "kcal/mol"

    assert info["N(QM-ATOMS)"] == 5
    assert units["N(QM-ATOMS)"] == "-"

    assert info["E(KIN)"] == 6
    assert units["E(KIN)"] == "kcal/mol"

    assert info["E(INTRA)"] == 7
    assert units["E(INTRA)"] == "kcal/mol"

    assert info["VOLUME"] == 8
    assert units["VOLUME"] == "A^3"

    assert info["DENSITY"] == 9
    assert units["DENSITY"] == "g/cm^3"

    assert info["MOMENTUM"] == 10
    assert units["MOMENTUM"] == "amuA/fs"

    assert info["LOOPTIME"] == 11
    assert units["LOOPTIME"] == "s"
