import pytest

from PQAnalysis.io import InfoFileReader
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.traj.exceptions import MDEngineFormatError
from PQAnalysis.exceptions import PQFileNotFoundError

from . import pytestmark



@pytest.mark.parametrize("example_dir", ["readInfoFile"], indirect=False)
def test__init__(test_with_data_dir):
    with pytest.raises(PQFileNotFoundError) as exception:
        InfoFileReader("tmp")
    assert str(exception.value) == "File tmp not found."

    with pytest.raises(MDEngineFormatError) as exception:
        InfoFileReader("md-01.info", engine_format="tmp")
    assert str(exception.value) == (
        "\n"
        "'tmp' is not a valid MDEngineFormat.\n"
        f"Possible values are: {MDEngineFormat.member_repr()} "
        "or their case insensitive string representation: "
        f"{MDEngineFormat.value_repr()}"
    )

    reader = InfoFileReader("md-01.info")
    assert reader.filename == "md-01.info"
    assert reader.format == MDEngineFormat.PQ

    reader = InfoFileReader("md-01.info", engine_format="qmcfc")
    assert reader.filename == "md-01.info"
    assert reader.format == MDEngineFormat.QMCFC



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

    reader = InfoFileReader("md-01.qmcfc.info", engine_format="qmcfc")
    info, units = reader.read()

    assert info["SIMULATION TIME"] == 0
    assert info["QM_MOLECULES"] == 1
    assert info["TEMPERATURE"] == 2
    assert info["PRESSURE"] == 3
    assert info["E(QM)"] == 4
    assert info["E(MM)"] == 5
    assert info["E(KIN)"] == 6
    assert info["E(INTRA)"] == 7
    assert info["E(BOND)"] == 8
    assert info["E(ANGLE)"] == 9
    assert info["E(DIHEDRAL)"] == 10
    assert info["E(IMPROPER)"] == 11
    assert info["E(COULOMB)"] == 12
    assert info["E(NONCOULOMB)"] == 13
    assert info["E(CF)"] == 14
    assert info["E(CF_RF)"] == 15
    assert info["E(RF)"] == 16
    assert info["E(THREEBODY)"] == 17
    assert info["E(NH_MOM)"] == 18
    assert info["E(NH_FRIC)"] == 19
    assert info["VOLUME"] == 20
    assert info["DENSITY"] == 21
    assert info["MOMENTUM"] == 22
    assert info["LOOPTIME"] == 23

    reader = InfoFileReader("md-01.qmcfc.info")
    with pytest.raises(MDEngineFormatError) as exception:
        reader.read()
    assert str(
        exception.value
    ) == "Info file md-01.qmcfc.info is not in PQ format."

    reader = InfoFileReader("md-01.info", engine_format="qmcfc")
    with pytest.raises(MDEngineFormatError) as exception:
        reader.read()
    assert str(
        exception.value
    ) == "Info file md-01.info is not in qmcfc format."
