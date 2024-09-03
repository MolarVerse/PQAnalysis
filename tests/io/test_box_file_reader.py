import pytest
import numpy as np

from . import pytestmark

from _pytest.capture import CaptureFixture

from PQAnalysis.io import BoxFileReader
from PQAnalysis.core.cell import Cell, Cells
from PQAnalysis.io import BoxWriter, write_box, BoxFileFormat, FileWritingMode
from PQAnalysis.io.exceptions import BoxWriterError, BoxFileFormatError, BoxReaderError
from PQAnalysis.traj import Trajectory
from PQAnalysis.core import Cell, Atom
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.exceptions import PQFileNotFoundError


class TestBoxFileReader:
    @pytest.mark.parametrize("example_dir", ["readBoxFile"], indirect=False)
    def test__init__(self, test_with_data_dir):
        with pytest.raises(PQFileNotFoundError) as exception:
            BoxFileReader(engine_format=MDEngineFormat.PQ)
        assert str(exception.value) == (
            "Either a filename or a trajectory must be "
            "provided depending on the engine format."
        )

        with pytest.raises(PQFileNotFoundError) as exception:
            BoxFileReader(filename="tmp", engine_format=MDEngineFormat.PQ)
        assert str(exception.value) == "File tmp not found."

        reader = BoxFileReader(filename="md-01.box",
                               engine_format=MDEngineFormat.PQ)
        assert reader.filename == "md-01.box"
        assert reader.engine_format == MDEngineFormat.PQ
        assert reader.trajectory is None

        trajectory = Trajectory()
        reader = BoxFileReader(trajectory=trajectory,
                               engine_format=MDEngineFormat.PQ)
        assert reader.trajectory == trajectory
        assert reader.engine_format == MDEngineFormat.PQ
        assert reader.filename is None

    @pytest.mark.parametrize("example_dir", ["readBoxFile"], indirect=False)
    def test__read_from_file(self, test_with_data_dir):
        reader = BoxFileReader(
            filename="md-01.box", engine_format=MDEngineFormat.PQ)
        cells = reader._read_from_file()
        assert cells == Cells([Cell(10.0, 10.0, 10.0), Cell(20.0, 20.0, 20.0)])

        reader = BoxFileReader(
            filename="md-02.box", engine_format=MDEngineFormat.PQ)
        cells = reader._read_from_file()
        assert cells == Cells([Cell(10.0, 10.0, 10.0, 90.0, 90.0, 120.0), Cell(
            20.0, 20.0, 20.0, 90.0, 90.0, 120.0)])

        reader = BoxFileReader(
            filename="md-invalid-columns.box", engine_format=MDEngineFormat.PQ)
        with pytest.raises(BoxReaderError) as exception:
            reader._read_from_file()
        assert str(exception.value) == (
            "Line 1 in file md-invalid-columns.box has "
            "an invalid number of columns: 5"
        )

        reader = BoxFileReader(
            filename="md-03.box", engine_format=MDEngineFormat.PQ)
        cells = reader._read_from_file()
        assert cells == Cells([Cell(10.0, 10.0, 10.0), Cell(
            20.0, 20.0, 20.0)])
