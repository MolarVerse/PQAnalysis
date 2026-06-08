import pytest
import numpy as np

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.analysis.rdf import RDFDataWriter, RDFLogWriter, RDF
from PQAnalysis.core import Atom, Cell
from PQAnalysis.traj import Trajectory
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

# import topology marker
from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception



class TestRDFLogWriter:

    def test__type_checking(self, caplog):
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
            function=RDFLogWriter,
            filename=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("rdf",
            1.0,
            RDF),
            exception=PQTypeError,
            function=RDFLogWriter("test.out").write_before_run,
            rdf=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("rdf",
            1.0,
            RDF),
            exception=PQTypeError,
            function=RDFLogWriter("test.out").write_after_run,
            rdf=1.0,
        )

    def test_write_before_run_writes_target_selection(self, tmp_path):
        traj = Trajectory([
            AtomicSystem(
                atoms=[Atom("H"), Atom("Q")],
                pos=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
                cell=Cell(10.0, 10.0, 10.0),
            )
        ])
        rdf = RDF(traj, "H", "Q", delta_r=1.0)
        log_file = tmp_path / "rdf.log"
        writer = RDFLogWriter(str(log_file))

        writer.write_before_run(rdf)
        writer.close()

        contents = log_file.read_text(encoding="utf-8")
        assert "    Target selection:    Q\n" in contents
        assert "{Selection" not in contents
