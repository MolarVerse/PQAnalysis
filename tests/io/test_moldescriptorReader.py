from pathlib import Path

import pytest
import numpy as np

from . import pytestmark

from PQAnalysis.io import MoldescriptorReader
from PQAnalysis.io.exceptions import MoldescriptorReaderError
from PQAnalysis.core import Element
from PQAnalysis.exceptions import PQFileNotFoundError



class TestMoldescriptorReader:

    @pytest.mark.parametrize(
        "example_dir",
        ["readMoldescriptor"],
        indirect=False
    )
    def test__init__(self, test_with_data_dir):
        reader = MoldescriptorReader("moldescriptor.dat")
        assert reader.filename == "moldescriptor.dat"

        with pytest.raises(PQFileNotFoundError) as exception:
            MoldescriptorReader("tmp")
        assert str(exception.value) == "File tmp not found."

    def test_read_mol_type(self):
        mol_type = f"""H2o 3 1.0
H 1 2.3
H 1 2.3
O 2 -4.6
"""

        mol_type = MoldescriptorReader._read_mol_type(mol_type.splitlines(), 1)

        assert mol_type.name == "H2o"
        assert mol_type.n_atoms == 3
        assert mol_type.total_charge == 1.0
        assert mol_type.elements == [Element("H"), Element("H"), Element("O")]
        assert np.allclose(mol_type.atom_types, np.array([1, 1, 2]))
        assert np.allclose(
            mol_type.partial_charges,
            np.array([2.3,
            2.3,
            -4.6])
        )
        assert mol_type.id == 1

        mol_type = f"""H2o 3 1.0
H 1 2.3
H 1
O 2 -4.6
"""

        with pytest.raises(MoldescriptorReaderError) as exception:
            MoldescriptorReader._read_mol_type(mol_type.splitlines(), 1)
        assert str(
            exception.value
        ) == "The number of columns in the body of a mol type must be 3 or 4.\nGot 2 columns instead in text: 'H 1'\n"

    @pytest.mark.parametrize(
        "example_dir",
        ["readMoldescriptor"],
        indirect=False
    )
    def test_read(self, test_with_data_dir):
        reader = MoldescriptorReader("moldescriptor.dat")
        mol_types = reader.read()

        assert len(mol_types) == 2
        water_type = mol_types[0]

        assert water_type.name == "H2O"
        assert water_type.n_atoms == 3
        assert np.allclose(water_type.total_charge, 0.0)
        assert water_type.elements == [
            Element("O"),
            Element("H"),
            Element("H")
        ]
        assert np.allclose(water_type.atom_types, np.array([0, 1, 1]))
        assert np.allclose(
            water_type.partial_charges,
            np.array([-0.65966,
            0.32983,
            0.32983])
        )
        assert water_type.id == 1

        ammonia_type = mol_types[1]

        assert ammonia_type.name == "Ammonia"
        assert ammonia_type.n_atoms == 4
        assert np.allclose(ammonia_type.total_charge, 0.0)
        assert ammonia_type.elements == [
            Element("N"),
            Element("H"),
            Element("H"),
            Element("H")
        ]
        assert np.allclose(ammonia_type.atom_types, np.array([0, 1, 1, 1]))
        assert np.allclose(
            ammonia_type.partial_charges,
            np.array([-0.8022,
            0.2674,
            0.2674,
            0.2674])
        )
        assert ammonia_type.id == 2

        reader = MoldescriptorReader("moldescriptor_withError.dat")
        with pytest.raises(MoldescriptorReaderError) as exception:
            reader.read()
        assert str(
            exception.value
        ) == "The number of columns in the header of a mol type must be 3.\nGot 2 columns instead in text: '  H2O            3'\n"

    @pytest.mark.parametrize(
        "example_dir",
        ["readMoldescriptor"],
        indirect=False
    )
    def test_read_too_few_atom_lines(self, test_with_data_dir):
        reader = MoldescriptorReader("moldescriptor_tooFewAtoms.dat")
        with pytest.raises(MoldescriptorReaderError) as exception:
            reader.read()
        assert str(exception.value) == (
            "The mol type H2O declares 3 atoms, but only 2 atom lines "
            "follow before the end of the file.\n"
        )

        reader = MoldescriptorReader("moldescriptor_tooFewAtoms_midFile.dat")
        with pytest.raises(MoldescriptorReaderError) as exception:
            reader.read()
        assert str(exception.value) == (
            "The mol type Ar declares 2 atoms, but only 0 atom lines "
            "follow before the end of the file.\n"
        )


WATER = "H2O 3 0.0\nO 0 -0.8\nH 1 0.4\nH 1 0.4\n"
AMMONIA = "NH3 4 0.0\nN 0 -0.8\nH 1 0.27\nH 1 0.27\nH 1 0.27\n"
ION = "Na 1 1.0\nNa 0 1.0\n"


def read(tmp_path, text):
    path = tmp_path / "moldescriptor.dat"
    path.write_text(text, encoding="utf-8")
    reader = MoldescriptorReader(str(path))

    return reader, reader.read()


class TestSpecialTypes:

    def test_the_types_are_kept(self, tmp_path):
        reader, mol_types = read(
            tmp_path,
            "WATER_TYPE 2\nammonia_type 3\n" + ION + WATER + AMMONIA,
        )

        assert (reader.water_type, reader.ammonia_type) == (2, 3)
        assert [m.name for m in mol_types] == ["Na", "H2O", "NH3"]
        assert mol_types[reader.water_type - 1].name == "H2O"
        assert mol_types[reader.ammonia_type - 1].name == "NH3"

    def test_without_the_lines_they_are_none(self, tmp_path):
        reader, mol_types = read(tmp_path, WATER)

        assert reader.water_type is None
        assert reader.ammonia_type is None
        assert len(mol_types) == 1

    def test_the_lines_can_follow_the_mol_types_and_have_comments(self, tmp_path):
        reader, _ = read(tmp_path, WATER + "WATER_TYPE 1 # the water\n")

        assert reader.water_type == 1

    def test_reading_again_resets_the_types(self, tmp_path):
        reader, _ = read(tmp_path, "WATER_TYPE 1\n" + WATER)
        Path(reader.filename).write_text(WATER, encoding="utf-8")

        reader.read()

        assert reader.water_type is None

    def test_the_returned_mol_types_do_not_depend_on_the_lines(self, tmp_path):
        _, with_lines = read(tmp_path, "WATER_TYPE 1\nAMMONIA_TYPE 2\n" + WATER + AMMONIA)
        _, without = read(tmp_path, WATER + AMMONIA)

        assert with_lines == without

    @pytest.mark.parametrize(
        "line", ["WATER_TYPE", "WATER_TYPE one", "WATER_TYPE 1 2", "AMMONIA_TYPE 1.5"]
    )
    def test_a_line_without_one_integer_raises(self, tmp_path, line):
        with pytest.raises(MoldescriptorReaderError, match="needs exactly one integer"):
            read(tmp_path, line + "\n" + WATER)

    def test_a_type_that_is_not_in_the_file_warns(self, tmp_path, caplog):
        with caplog.at_level("WARNING"):
            reader, _ = read(tmp_path, "WATER_TYPE 3\n" + WATER)

        assert reader.water_type == 3
        assert "WATER_TYPE 3 points to a mol type that is not in the file" in " ".join(
            caplog.text.split()
        )

    def test_ammonia_type_not_in_the_file_warns(self, tmp_path, caplog):
        with caplog.at_level("WARNING"):
            read(tmp_path, "AMMONIA_TYPE 2\n" + WATER)

        assert "AMMONIA_TYPE 2 points to a mol type that is not in the file" in " ".join(
            caplog.text.split()
        )

    def test_a_water_type_that_is_no_water_warns(self, tmp_path, caplog):
        with caplog.at_level("WARNING"):
            read(tmp_path, "WATER_TYPE 1\n" + AMMONIA)

        assert "not made of the atoms O, H and H in this order" in " ".join(
            caplog.text.split()
        )

    def test_a_correct_file_logs_no_warning(self, tmp_path, caplog):
        with caplog.at_level("WARNING"):
            read(tmp_path, "WATER_TYPE 1\nAMMONIA_TYPE 2\n" + WATER + AMMONIA)

        assert caplog.text == ""

    @pytest.mark.parametrize("example_dir", ["readMoldescriptor"], indirect=False)
    def test_the_example_file(self, test_with_data_dir):
        reader = MoldescriptorReader("moldescriptor.dat")
        reader.read()

        assert (reader.water_type, reader.ammonia_type) == (1, 2)
