import pytest
import numpy as np

from . import pytestmark

from PQAnalysis.io.virial import VirialFileReader, StressFileReader
from PQAnalysis.io.virial.api import read_virial_file, read_stress_file
from PQAnalysis.io.exceptions import VirialFileReaderError



class TestVirialFileReader:

    @pytest.mark.usefixtures("tmpdir")
    def test_read_skips_blank_and_comment_lines(self):
        with open("md.vir", "w", encoding="utf-8") as file:
            print("# step xx xy xz yx yy yz zx zy zz", file=file)
            print("1 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0", file=file)
            print("", file=file)

        data = VirialFileReader("md.vir").read()

        assert len(data) == 1
        assert np.allclose(
            data[0],
            np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        )

        data = read_virial_file("md.vir")
        assert len(data) == 1

    @pytest.mark.usefixtures("tmpdir")
    def test_read_invalid_number_of_columns(self):
        with open("md.vir", "w", encoding="utf-8") as file:
            print("1 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0", file=file)
            print("2 1.0 2.0 3.0", file=file)

        with pytest.raises(VirialFileReaderError) as exception:
            VirialFileReader("md.vir").read()
        assert str(exception.value) == (
            "Invalid number of columns in file md.vir line 2. "
            "Expected 10 columns."
        )



class TestStressFileReader:

    @pytest.mark.usefixtures("tmpdir")
    def test_read_skips_blank_and_comment_lines(self):
        with open("md.stress", "w", encoding="utf-8") as file:
            print("# step xx xy xz yx yy yz zx zy zz", file=file)
            print("1 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0", file=file)
            print("", file=file)

        data = StressFileReader("md.stress").read()

        assert len(data) == 1
        assert np.allclose(
            data[0],
            np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        )

        data = read_stress_file("md.stress")
        assert len(data) == 1

    @pytest.mark.usefixtures("tmpdir")
    def test_read_invalid_number_of_columns(self):
        with open("md.stress", "w", encoding="utf-8") as file:
            print("1 1.0 2.0 3.0", file=file)

        with pytest.raises(VirialFileReaderError) as exception:
            StressFileReader("md.stress").read()
        assert str(exception.value) == (
            "Invalid number of columns in file md.stress line 1. "
            "Expected 10 columns."
        )
