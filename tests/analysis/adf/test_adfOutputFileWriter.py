"""
A module to test the ADF output/log file writers.
"""

import numpy as np
import pytest  # pylint: disable=unused-import

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.analysis.adf import ADFDataWriter, ADFLogWriter, ADF
from PQAnalysis.core import Atom, Cell
from PQAnalysis.traj import Trajectory
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

# import topology marker
from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception



def _make_adf():
    traj = Trajectory([
        AtomicSystem(
            atoms=[Atom("O"), Atom("H"), Atom("H")],
            pos=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            cell=Cell(10.0, 10.0, 10.0),
        )
    ])
    return ADF(traj, "O", "H", n_angle_bins=180)



class TestADFDataWriter:

    def test_type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("filename", 1.0, str),
            exception=PQTypeError,
            function=ADFDataWriter,
            filename=1.0,
        )

    def test_write_four_columns(self, tmp_path):
        adf = _make_adf()
        data = adf.run()

        out_file = tmp_path / "adf.out"
        writer = ADFDataWriter(str(out_file))
        writer.write(data)

        lines = [
            line for line in out_file.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]

        assert len(lines) == adf.n_angle_bins
        assert all(len(line.split()) == 4 for line in lines)



class TestADFLogWriter:

    def test_type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "filename",
            1.0,
            str | None
            ),
            exception=PQTypeError,
            function=ADFLogWriter,
            filename=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("adf", 1.0, ADF),
            exception=PQTypeError,
            function=ADFLogWriter("test.out").write_before_run,
            adf=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("adf", 1.0, ADF),
            exception=PQTypeError,
            function=ADFLogWriter("test.out").write_after_run,
            adf=1.0,
        )

    def test_write_before_run_writes_selections(self, tmp_path):
        adf = _make_adf()
        log_file = tmp_path / "adf.log"
        writer = ADFLogWriter(str(log_file))

        writer.write_before_run(adf)
        writer.close()

        contents = log_file.read_text(encoding="utf-8")
        assert "ADF calculation:" in contents
        assert "Number of angle bins: 180" in contents
        assert "Reference selection: O" in contents
        assert "Target selection (ligand 1): H" in contents
        assert "Target selection (ligand 2): H" in contents
        assert "i-j radial gate: none" in contents
        assert "{Selection" not in contents

    def test_write_before_run_writes_gates(self, tmp_path):
        traj = Trajectory([
            AtomicSystem(
                atoms=[Atom("O"), Atom("H"), Atom("H")],
                pos=np.array(
                    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
                ),
                cell=Cell(10.0, 10.0, 10.0),
            )
        ])
        adf = ADF(traj, "O", "H", n_angle_bins=90, r_min_1=0.8, r_max_1=1.2)

        log_file = tmp_path / "adf.log"
        writer = ADFLogWriter(str(log_file))
        writer.write_before_run(adf)
        writer.close()

        contents = log_file.read_text(encoding="utf-8")
        assert "i-j radial gate: [0.8, 1.2)" in contents
        assert "i-k radial gate: none" in contents

    def test_write_after_run_writes_elapsed_time(self, tmp_path):
        adf = _make_adf()
        adf.run()

        log_file = tmp_path / "adf.log"
        writer = ADFLogWriter(str(log_file))
        writer.write_after_run(adf)
        writer.close()

        assert "Elapsed time:" in log_file.read_text(encoding="utf-8")
