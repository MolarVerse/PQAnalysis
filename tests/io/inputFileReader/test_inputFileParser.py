import pytest

from PQAnalysis.io import InputFileParser


class TestInputFileParser:
    @pytest.mark.parametrize("example_dir", ["inputFileReader"], indirect=False)
    def test__init__(self, test_with_data_dir):
        with pytest.raises(FileNotFoundError) as exception:
            InputFileParser("non-existent-file")
        assert str(
            exception.value) == "File non-existent-file not found."
