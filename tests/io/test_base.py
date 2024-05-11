import pytest
import os
import sys

from . import pytestmark

from PQAnalysis.io import BaseWriter, BaseReader, FileWritingMode
from PQAnalysis.io.exceptions import FileWritingModeError
from PQAnalysis.exceptions import PQFileNotFoundError



class TestBaseWriter:

    @pytest.mark.usefixtures("tmpdir")
    def test__init__(self):

        filename = "tmp"

        with pytest.raises(FileWritingModeError) as exception:
            BaseWriter(filename, "r")
            assert str(
                exception.value
            ) == """"
'r' is not a valid FileWritingMode.
Possible values are: FileWritingMode.OVERWRITE, FileWritingMode.APPEND, FileWritingMode.WRITE
or their case insensitive string representation: o, a, w"""

        open(filename, "w")

        with pytest.raises(FileWritingModeError) as exception:
            BaseWriter(filename, "w")
        assert str(
            exception.value
        ) == "File tmp already exists. Use mode \'a\' to append to the file or mode 'o' to overwrite the file."

        writer = BaseWriter(filename, "o")
        assert writer.file is None
        assert writer.mode == FileWritingMode.WRITE
        assert writer.filename == filename

        writer = BaseWriter(filename, "a")
        assert writer.file is None
        assert writer.mode == FileWritingMode.APPEND
        assert writer.filename == filename

        os.remove(filename)

        writer = BaseWriter(filename, "w")
        assert writer.file is None
        assert writer.mode == FileWritingMode.WRITE
        assert writer.filename == filename

        writer = BaseWriter(filename, "o")
        assert writer.file is None
        assert writer.mode == FileWritingMode.WRITE
        assert writer.filename == filename

        writer = BaseWriter()
        assert writer.file == sys.stdout
        assert writer.mode == FileWritingMode.WRITE
        assert writer.filename is None

    @pytest.mark.usefixtures("tmpdir")
    def test_open(self):

        filename = "tmp"

        writer = BaseWriter(filename, "a")
        assert writer.file is None

        writer.open()
        assert writer.file is not None

        os.remove(filename)

        writer = BaseWriter(filename, "w")
        assert writer.file is None

        writer.open()
        assert writer.file is not None

    @pytest.mark.usefixtures("tmpdir")
    def test_close(self):

        filename = "tmp"

        writer = BaseWriter(filename, "a")
        assert writer.file is None

        writer.open()
        assert writer.file is not None

        writer.close()
        assert writer.file is None



class TestBaseReader:

    @pytest.mark.usefixtures("tmpdir")
    def test__init__(self):
        with pytest.raises(PQFileNotFoundError) as exception:
            BaseReader("tmp")
        assert str(exception.value) == "File tmp not found."

        filename = "tmp"
        file = open(filename, "w")
        print("test", file=file)
        file.close()

        reader = BaseReader(filename)
        assert reader.filename == filename
        assert reader.multiple_files is False

        with pytest.raises(PQFileNotFoundError) as exception:
            BaseReader([filename, "tmp2"])
        assert str(
            exception.value
        ) == "At least one of the given files does not exist. File tmp2 not found."

        filename2 = "tmp2"
        file = open(filename2, "w")
        print("test2", file=file)
        file.close()

        reader = BaseReader([filename, filename2])
        assert reader.filenames == [filename, filename2]
        assert reader.multiple_files is True
