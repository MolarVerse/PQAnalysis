import logging
from pathlib import Path

import pytest

from PQAnalysis.analysis.rdf.rdf_input_file_reader import RDFInputFileReader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError, InputFileWarning
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError, PQFileNotFoundError

# import topology marker
from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception



class TestRDFInputFileReader:

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test__init__(self, test_with_data_dir):
        with pytest.raises(PQFileNotFoundError) as exception:
            RDFInputFileReader("not-a-file")
        assert str(exception.value) == "File not-a-file not found."

        filename = "input.in"

        reader = RDFInputFileReader(filename)
        assert reader.filename == filename
        assert reader.parser.filename == filename

    def test__init__type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("filename", 1.0, str),
            exception=PQTypeError,
            function=RDFInputFileReader,
            filename=1.0,
        )

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_restart_no_intra_molecular_not_none(self, test_with_data_dir):
        reader = RDFInputFileReader("input_no_intra_molecular_restart.in")
        with pytest.raises(InputFileError) as exception:
            reader.read()
        assert str(exception.value) == (
            "The no_intra_molecular key requires both a restart file and a "
            "moldescriptor file. Could not infer missing files from trajectory "
            "file 'md-00.xyz'. Tried restart_file='md-00.rst' and "
            "moldescriptor_file='moldescriptor.dat'."
        )

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_no_intra_molecular_false_does_not_require_topology(
        self, test_with_data_dir
    ):
        reader = RDFInputFileReader("input_no_intra_molecular_false.in")

        reader.read()

        assert reader.no_intra_molecular is False

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    @pytest.mark.parametrize(
        "keys, expected",
        [
            ("reference_com = True\ntarget_com = True\n", (True, True)),
            ("reference_com = False\n", (False, None)),
            ("target_com = True\n", (None, True)),
            ("", (None, None)),
        ],
    )
    def test_center_of_mass_keys(self, test_with_data_dir, keys, expected):
        base = Path("input_topology_files_default.in").read_text(
            encoding="utf-8"
        )
        Path("input_com.in").write_text(base + keys, encoding="utf-8")
        reader = RDFInputFileReader("input_com.in")

        reader.read()

        assert (reader.reference_com, reader.target_com) == expected

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    @pytest.mark.parametrize(
        "keys", ["reference_com = True", "target_com = True"]
    )
    def test_center_of_mass_keys_need_the_topology_files(
        self, test_with_data_dir, keys
    ):
        base = Path("input_no_restart.in").read_text(encoding="utf-8")
        Path("input_com.in").write_text(base + "\n" + keys + "\n", encoding="utf-8")
        reader = RDFInputFileReader("input_com.in")

        with pytest.raises(InputFileError) as exception:
            reader.read()

        assert str(exception.value).startswith(
            "The reference_com or target_com key requires both a restart "
            "file and a moldescriptor file. Could not infer missing files "
            "from trajectory file"
        )

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_center_of_mass_keys_infer_the_topology_files(
        self, test_with_data_dir
    ):
        Path("moldescriptor.dat").write_text(
            "C 1 0.0\nC 0 0.0\n", encoding="utf-8"
        )
        base = Path("input_no_intra_molecular_infer.in").read_text(
            encoding="utf-8"
        ).replace("no_intra_molecular = True", "")
        Path("input_com.in").write_text(
            base + "\nreference_com = True\ntarget_com = True\n",
            encoding="utf-8",
        )
        reader = RDFInputFileReader("input_com.in")

        reader.read()

        assert reader.restart_file == "infer.rst"
        assert reader.moldescriptor_file == "moldescriptor.dat"
        assert reader.reference_com is True

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_center_of_mass_false_does_not_require_topology(
        self, test_with_data_dir
    ):
        base = Path("input_no_restart.in").read_text(encoding="utf-8")
        Path("input_com.in").write_text(
            base + "\nreference_com = False\ntarget_com = False\n",
            encoding="utf-8",
        )

        RDFInputFileReader("input_com.in").read()

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_moldescriptor_no_restart(self, test_with_data_dir):
        reader = RDFInputFileReader("input_moldescriptor_restart.in")
        with pytest.raises(InputFileError) as exception:
            reader.read()
        assert str(exception.value) == (
            "The moldescriptor_file key can only be used in a meaningful way "
            "if a restart file is given."
        )

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_no_intra_molecular_infers_topology_files(
        self, test_with_data_dir, capture_log
    ):
        Path("moldescriptor.dat").write_text(
            "C 1 0.0\nC 0 0.0\n",
            encoding="utf-8",
        )
        reader = RDFInputFileReader("input_no_intra_molecular_infer.in")

        with capture_log(logging.INFO, RDFInputFileReader.logger) as log:
            reader.read()

        assert reader.no_intra_molecular is True
        assert reader.restart_file == "infer.rst"
        assert reader.moldescriptor_file == "moldescriptor.dat"
        assert "Inferred restart_file='infer.rst' from traj_files." in log.text
        assert (
            "Inferred moldescriptor_file='moldescriptor.dat' from traj_files."
            in log.text
        )

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_topology_files_default_to_no_intra_molecular(
        self, test_with_data_dir, capture_log
    ):
        reader = RDFInputFileReader("input_topology_files_default.in")

        with capture_log(logging.INFO, RDFInputFileReader.logger) as log:
            reader.read()

        assert reader.no_intra_molecular is True
        assert "Enabled no_intra_molecular because both restart_file" in log.text
        assert "moldescriptor_file are set." in log.text

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_explicit_no_intra_molecular_false_is_preserved(
        self, test_with_data_dir
    ):
        reader = RDFInputFileReader("input_topology_files_false.in")

        reader.read()

        assert reader.no_intra_molecular is False

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_no_intra_molecular_empty_traj_glob_reports_input_error(
        self, test_with_data_dir
    ):
        reader = RDFInputFileReader("input_no_intra_molecular_empty_glob.in")

        with pytest.raises(InputFileError) as exception:
            reader.read()

        assert str(exception.value) == (
            "The no_intra_molecular key requires both a restart file and a "
            "moldescriptor file. Could not infer missing files because no "
            "trajectory files were available."
        )



class TestRDFInputFileReaderPathResolution:

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_inference_from_subdirectory(self, test_with_data_dir):
        """
        Companion files are inferred beside the trajectory as written in
        the input file and resolved once, not twice.
        """
        run = Path("run")
        run.mkdir()
        for name in ("input_no_intra_molecular_infer.in", "infer.xyz", "infer.rst"):
            Path(name).rename(run / name)
        (run / "moldescriptor.dat").write_text("C 1 0.0\nC 0 0.0\n", encoding="utf-8")

        reader = RDFInputFileReader("run/input_no_intra_molecular_infer.in")
        reader.read()

        assert reader.no_intra_molecular is True
        assert reader.traj_files == ["run/infer.xyz"]
        assert reader.restart_file == "run/infer.rst"
        assert reader.moldescriptor_file == "run/moldescriptor.dat"
        assert reader.out_file == "run/rdf.out"
