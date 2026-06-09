import pytest
import numpy as np

from . import pytestmark

from PQAnalysis.exceptions import PQFileNotFoundError
from PQAnalysis.io import BoxReader, read_box
from PQAnalysis.io.exceptions import BoxReaderError



class TestBoxReader:

    @pytest.mark.usefixtures("tmpdir")
    def test__init__(self):
        with pytest.raises(PQFileNotFoundError) as exception:
            BoxReader("tmp")
        assert str(exception.value) == "File tmp not found."

        with open("box.dat", "w", encoding="utf-8") as file:
            print("1 1.0 2.0 3.0 90.0 90.0 90.0", file=file)

        reader = BoxReader("box.dat")

        assert reader.filename == "box.dat"
        assert reader.multiple_files is False

    @pytest.mark.usefixtures("tmpdir")
    def test_read(self):
        with open("box.dat", "w", encoding="utf-8") as file:
            print("# step x y z alpha beta gamma", file=file)
            print("1 1.1 1.2 1.3 90.0 91.0 92.0", file=file)
            print("", file=file)
            print("2 2.1 2.2 2.3 80.0 81.0 82.0", file=file)

        steps, box_lengths, box_angles = BoxReader("box.dat").read()

        assert np.array_equal(steps, np.array([1, 2]))
        assert np.allclose(
            box_lengths,
            np.array([[1.1, 1.2, 1.3], [2.1, 2.2, 2.3]])
        )
        assert np.allclose(
            box_angles,
            np.array([[90.0, 91.0, 92.0], [80.0, 81.0, 82.0]])
        )

    @pytest.mark.usefixtures("tmpdir")
    def test_read_box(self):
        with open("box.dat", "w", encoding="utf-8") as file:
            print("1 1.0 2.0 3.0 90.0 91.0 92.0", file=file)

        steps, box_lengths, box_angles = read_box("box.dat")

        assert np.array_equal(steps, np.array([1]))
        assert np.allclose(box_lengths, np.array([[1.0, 2.0, 3.0]]))
        assert np.allclose(box_angles, np.array([[90.0, 91.0, 92.0]]))

    @pytest.mark.usefixtures("tmpdir")
    def test_invalid_number_of_columns(self):
        with open("box.dat", "w", encoding="utf-8") as file:
            print("1 1.0 1.0", file=file)

        with pytest.raises(BoxReaderError) as exception:
            BoxReader("box.dat").read()
        assert str(exception.value) == (
            "Invalid number of columns in box file line 1. Expected 7 columns."
        )

    @pytest.mark.usefixtures("tmpdir")
    def test_invalid_numeric_value(self):
        with open("box.dat", "w", encoding="utf-8") as file:
            print("1 bad 1.0 1.0 90.0 90.0 90.0", file=file)

        with pytest.raises(BoxReaderError) as exception:
            BoxReader("box.dat").read()
        assert str(exception.value) == (
            "Invalid numeric value in box file line 1: "
            "1 bad 1.0 1.0 90.0 90.0 90.0"
        )
