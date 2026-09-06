import pytest

from PQAnalysis.io.base import PQFileNotFoundError

from ... import pytestmark

from PQAnalysis.io.input_file_reader.pq_analysis.pqanalysis_input_file_reader import PQAnalysisInputFileReader
from PQAnalysis.io.input_file_reader.formats import InputFileFormat

from PQAnalysis.io.input_file_reader.exceptions import InputFileError



class TestPQAnalysisInputFileReader:

    @pytest.mark.parametrize(
        "example_dir", ["inputFileReader/PQAnalysis_input/"], indirect=False
    )
    def test__init__(self, test_with_data_dir):
        pqanalysis_input_file_reader = PQAnalysisInputFileReader("input_PQ.in")

        assert pqanalysis_input_file_reader.filename == "input_PQ.in"
        assert pqanalysis_input_file_reader.format == InputFileFormat.PQANALYSIS
        assert pqanalysis_input_file_reader.parser.filename == "input_PQ.in"
        assert pqanalysis_input_file_reader.parser.input_format == InputFileFormat.PQANALYSIS

        assert pqanalysis_input_file_reader.dictionary is None
        assert pqanalysis_input_file_reader.raw_input_file is None

    def test__init__no_input_file(self):
        with pytest.raises(PQFileNotFoundError) as exception:
            PQAnalysisInputFileReader("no_input_file.in")
        assert str(exception.value) == f"File no_input_file.in not found."

    @pytest.mark.parametrize(
        "example_dir", ["inputFileReader/PQAnalysis_input/"], indirect=False
    )
    def test_read(self, test_with_data_dir):
        pqanalysis_input_file_reader = PQAnalysisInputFileReader("input_PQ.in")

        pqanalysis_input_file_reader.read()

        assert pqanalysis_input_file_reader.raw_input_file == open(
            "input_PQ.in", "r"
        ).read()
        assert pqanalysis_input_file_reader.dictionary == pqanalysis_input_file_reader.parser.parse(
        )

    @pytest.mark.parametrize(
        "example_dir", ["inputFileReader/PQAnalysis_input/"], indirect=False
    )
    def test_check_required_keys(self, test_with_data_dir):
        pqanalysis_input_file_reader = PQAnalysisInputFileReader("input_PQ.in")
        pqanalysis_input_file_reader.read()

        required_keys = ["required_key1", "required_key2"]

        assert pqanalysis_input_file_reader.check_required_keys(
            required_keys
        ) == None

        required_keys = ["required_key1", "required_key2", "required_key3"]
        with pytest.raises(InputFileError) as exception:
            pqanalysis_input_file_reader.check_required_keys(required_keys)
        assert str(exception.value) == (
            "Not all required keys were set in the input file! The required keys are: ['required_key1', 'required_key2', 'required_key3']."
        )

    @pytest.mark.parametrize(
        "example_dir", ["inputFileReader/PQAnalysis_input/"], indirect=False
    )
    def test_check_known_keys(self, test_with_data_dir):
        pqanalysis_input_file_reader = PQAnalysisInputFileReader("input_PQ.in")
        pqanalysis_input_file_reader.read()

        known_keys = [
            "required_key1", "required_key2", "optional_key1", "optional_key2"
        ]

        assert pqanalysis_input_file_reader.check_known_keys(
            known_keys
        ) == None

    @pytest.mark.parametrize(
        "example_dir", ["inputFileReader/PQAnalysis_input/"], indirect=False
    )
    def test_check_known_keys_no_known_keys(self, test_with_data_dir):
        pqanalysis_input_file_reader = PQAnalysisInputFileReader("input_PQ.in")
        pqanalysis_input_file_reader.read()
        known_keys = None
        assert pqanalysis_input_file_reader.check_known_keys(
            known_keys
        ) == None

    @pytest.mark.parametrize(
        "example_dir", ["inputFileReader/PQAnalysis_input/"], indirect=False
    )
    def test_not_defined_optional_keys(self, test_with_data_dir):
        pqanalysis_input_file_reader = PQAnalysisInputFileReader("input_PQ.in")
        pqanalysis_input_file_reader.read()

        optional_keys = ["optional_key1", "optional_key2", "optional_key3"]

        assert pqanalysis_input_file_reader.not_defined_optional_keys(
            optional_keys
        ) == None
        assert pqanalysis_input_file_reader.dictionary["optional_key3"] == (
            None, "None", "None"
        )



class TestPathResolution:

    def _write(self, path, text):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def test_relative_paths_resolve_against_input_file_directory(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        self._write(
            tmp_path / "run" / "input.in",
            "traj_files = traj.xyz\n"
            "out_file = rdf.dat\n"
            "log_file = ../rdf.log\n"
            "restart_file = data/start.rst\n",
        )

        reader = PQAnalysisInputFileReader("run/input.in")
        reader.read()

        assert reader.traj_files == ["run/traj.xyz"]
        assert reader.out_file == "run/rdf.dat"
        assert reader.log_file == "run/../rdf.log"
        assert reader.restart_file == "run/data/start.rst"
        assert reader.moldescriptor_file is None

    def test_input_file_in_working_directory_leaves_paths_unchanged(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        self._write(
            tmp_path / "input.in",
            "traj_files = traj.xyz\nout_file = rdf.dat\n",
        )

        reader = PQAnalysisInputFileReader("input.in")
        reader.read()

        assert reader.traj_files == ["traj.xyz"]
        assert reader.out_file == "rdf.dat"

    def test_absolute_paths_are_not_changed(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        absolute = tmp_path / "elsewhere" / "traj.xyz"
        self._write(
            tmp_path / "run" / "input.in",
            f"traj_files = {absolute}\nout_file = {tmp_path / 'rdf.dat'}\n",
        )

        reader = PQAnalysisInputFileReader("run/input.in")
        reader.read()

        assert reader.traj_files == [str(absolute)]
        assert reader.out_file == str(tmp_path / "rdf.dat")

    def test_glob_expands_in_input_file_directory(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        for name in ("md-01.xyz", "md-02.xyz"):
            self._write(tmp_path / "run" / name, "")
        self._write(tmp_path / "md-99.xyz", "")  # in cwd, must not match
        self._write(tmp_path / "run" / "input.in", "traj_files = md-*.xyz\n")

        reader = PQAnalysisInputFileReader("run/input.in")
        reader.read()

        assert sorted(reader.traj_files) == ["run/md-01.xyz", "run/md-02.xyz"]

    def test_resolve_path_none(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        self._write(tmp_path / "input.in", "traj_files = traj.xyz\n")
        reader = PQAnalysisInputFileReader("input.in")

        assert reader.resolve_path(None) is None
        assert reader.resolve_paths(None) is None
