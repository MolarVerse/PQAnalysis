import pytest
import numpy as np

from . import pytestmark

from PQAnalysis.io import RestartFileReader
from PQAnalysis.io.restart_file.exceptions import RestartFileReaderError
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.core import Atom, Cell, Residue, Element
from PQAnalysis.exceptions import PQFileNotFoundError



class Test_RestartFileReader:

    @pytest.mark.parametrize(
        "example_dir", ["readRestartFile"], indirect=False
    )
    def test__init__(self, test_with_data_dir):
        with pytest.raises(PQFileNotFoundError) as exception:
            RestartFileReader("tmp")
        assert str(exception.value) == "File tmp not found."

        filename = "md-01.rst"

        reader = RestartFileReader(filename)
        assert reader.filename == filename
        assert reader.md_engine_format == MDEngineFormat.PQ
        assert reader.moldescriptor_filename == None
        assert reader.reference_residues == None

        reader = RestartFileReader(
            filename,
            md_engine_format="qmcfc",
            moldescriptor_filename="moldescriptor.dat",
        )
        assert reader.filename == filename
        assert reader.md_engine_format == MDEngineFormat.QMCFC
        assert reader.moldescriptor_filename == "moldescriptor.dat"
        assert reader.reference_residues == None

        residue = Residue("H2O", 1, 0, Element("Ar"), 1, 0)
        reader = RestartFileReader(
            filename, md_engine_format="qmcfc", reference_residues=[residue]
        )
        assert reader.filename == filename
        assert reader.md_engine_format == MDEngineFormat.QMCFC
        assert reader.moldescriptor_filename == None
        assert reader.reference_residues == [residue]

        with pytest.raises(RestartFileReaderError) as exception:
            RestartFileReader(
                filename,
                md_engine_format="qmcfc",
                moldescriptor_filename="moldescriptor.dat",
                reference_residues=[residue]
            )
        assert str(
            exception.value
        ) == "Both moldescriptor_filename and reference_residues are given. They are mutually exclusive."

    def test__parse_box(self):
        line = ["box"]
        assert RestartFileReader._parse_box(line) == Cell()

        line = ["box", "1.0", "2.0", "3.0"]
        assert RestartFileReader._parse_box(line) == Cell(1.0, 2.0, 3.0)

        line = ["box", "1.0", "2.0", "3.0", "90.0", "90.0", "120.0"]
        assert RestartFileReader._parse_box(line) == Cell(
            1.0, 2.0, 3.0, 90.0, 90.0, 120.0
        )

        with pytest.raises(RestartFileReaderError) as exception:
            line = ["box", "1.0", "2.0"]
            RestartFileReader._parse_box(line)
        assert str(exception.value) == "Invalid number of arguments for box: 3"

    def test__parse_atoms(self):
        with pytest.raises(RestartFileReaderError) as exception:
            lines = []
            RestartFileReader._parse_atoms(lines, Cell())
        assert str(exception.value) == "No atoms found in restart file."

        with pytest.raises(RestartFileReaderError) as exception:
            lines = ["C 0 1 1.0 1.0 1.0"]
            RestartFileReader._parse_atoms(lines, Cell())
        assert str(
            exception.value
        ) == "Invalid number of arguments for atom: 6"

        lines = [
            "C 0 1 1.0 1.0 1.0 1.1 1.2 1.3 1.4 1.5 1.6",
            "H 0 2 2.0 2.0 2.0 2.1 2.2 2.3 2.4 2.5 2.6"
        ]
        frame = RestartFileReader._parse_atoms(lines, Cell())
        system = frame
        residue_ids = frame.topology.residue_ids
        assert system.n_atoms == 2
        assert system.atoms == [Atom(name) for name in ["C", "H"]]
        assert system.pos.shape == (2, 3)
        assert np.allclose(
            system.pos, np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
        )
        assert system.vel.shape == (2, 3)
        assert np.allclose(
            system.vel, np.array([[1.1, 1.2, 1.3], [2.1, 2.2, 2.3]])
        )
        assert system.forces.shape == (2, 3)
        assert np.allclose(
            system.forces, np.array([[1.4, 1.5, 1.6], [2.4, 2.5, 2.6]])
        )
        assert np.allclose(residue_ids, np.array([1, 2]))
        assert system.cell == Cell()

    @pytest.mark.parametrize(
        "example_dir", ["readRestartFile"], indirect=False
    )
    def test_read(self, test_with_data_dir):
        frame = RestartFileReader(
            "md-01.rst", moldescriptor_filename="moldescriptor.dat"
        ).read()
        system = frame
        residue_ids = frame.topology.residue_ids

        assert np.allclose(residue_ids, np.array([1, 1, 2, 2]))
        assert system.n_atoms == 4
        assert system.atoms == [Atom(name) for name in ["C", "H", "N", "N"]]

        assert system.pos.shape == (4, 3)
        assert np.allclose(
            system.pos,
            np.array(
                [
                    [1.0, 1.1, 1.2], [2.0, 2.1, 2.2], [3.0, 3.1, 3.2],
                    [4.0, 4.1, 4.2]
                ]
            )
        )
        assert system.vel.shape == (4, 3)
        assert np.allclose(
            system.vel,
            np.array(
                [
                    [1.3, 1.4, 1.5], [2.3, 2.4, 2.5], [3.3, 3.4, 3.5],
                    [4.3, 4.4, 4.5]
                ]
            )
        )
        assert system.forces.shape == (4, 3)
        assert np.allclose(
            system.forces,
            np.array(
                [
                    [1.6, 1.7, 1.8], [2.6, 2.7, 2.8], [3.6, 3.7, 3.8],
                    [4.6, 4.7, 4.8]
                ]
            )
        )
        assert system.cell == Cell(
            15.0623, 15.0964, 20.0232, 89.9232, 90.2261, 120.324
        )
        assert np.allclose(system.topology.residue_ids, np.array([1, 1, 2, 2]))
        assert len(system.topology.reference_residues) == 2
