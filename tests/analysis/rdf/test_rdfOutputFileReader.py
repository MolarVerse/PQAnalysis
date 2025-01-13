import pytest
import numpy as np
from beartype.typing import Tuple

from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.io.traj_file.trajectory_reader import TrajectoryReader
from PQAnalysis.topology.topology import Topology
from PQAnalysis.analysis.rdf import RDFDataWriter, RDFLogWriter, RDF

from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError

# import topology marker
from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception


class TestRDFDataWriter:

    def test__type_checking(self,caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "filename",
            1.0,
            str
            ),
            exception=PQTypeError,
            function=RDFDataWriter,
            filename=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
            "data",
            1.0,
            Tuple[Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray]
            ),
            exception=PQTypeError,
            function=RDFDataWriter("test.out").write,
            data=1.0,
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "data",
                Tuple[np.array([1.0, 2.0, 3.0]),
                np.array([4.0, 5.0, 6.0]),
                np.array([7.0, 8.0, 9.0]),
                np.array([10.0, 11.0, 12.0]),
                np.array([13.0, 14.0, 15.0])],
                Tuple[Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray]
            ),
            exception=PQTypeError,
            function=RDFDataWriter("test.out").write,
            data=Tuple[np.array([1.0, 2.0, 3.0]),
                np.array([4.0, 5.0, 6.0]),
                np.array([7.0, 8.0, 9.0]),
                np.array([10.0, 11.0, 12.0]),
                np.array([13.0, 14.0, 15.0])],
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "data",
                Tuple[[1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0],
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]],
                Tuple[Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray]
            ),
            exception=PQTypeError,
            function=RDFDataWriter("test.out").write,
            data=tuple[[1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0],
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]],
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "data",
                [[1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0],
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]],
                Tuple[Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray]
            ),
            exception=PQTypeError,
            function=RDFDataWriter("test.out").write,
            data=[[1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0],
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]],
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message(
                "data",
                tuple(np.array([1.0, 2.0, 3.0])),
                Tuple[Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray,
        Np1DNumberArray]
            ),
            exception=PQTypeError,
            function=RDFDataWriter("test.out").write,
            data=tuple(np.array([1.0, 2.0, 3.0]))
        )

    def test__init__(self):
        writer = RDFDataWriter("test.out")
        assert writer.filename == "test.out"
        writer.close()
    
    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_write(self, test_with_data_dir):
        writer = RDFDataWriter("test.out")
        data = tuple([np.array([1.0, 2.0, 3.0]),np.array([4.0, 5.0, 6.0]),np.array([7.0, 8.0, 9.0]),np.array([10.0, 11.0, 12.0]), np.array([13.0, 14.0, 15.0])])
        writer.write(
            (
                data
            )
        )
        
        writer.close()

        with open("test.out") as file:
            lines = file.readlines()

        assert lines == [
            "1.0 4.0 7.0 10.0 13.0\n",
            "2.0 5.0 8.0 11.0 14.0\n",
            "3.0 6.0 9.0 12.0 15.0\n",
        ]

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

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_write_before_run(self, test_with_data_dir):
        writer = RDFLogWriter("test.out")
        traj_reader = TrajectoryReader("traj.xyz",topology=Topology())
        ref_selection = "H"
        target_selection = "O"
        rdf = RDF(
            traj=traj_reader,
            reference_species=ref_selection,
            target_species=target_selection,
            delta_r=0.1
        )

        writer.write_before_run(rdf)
        writer.close()

        with open("test.out",encoding="utf-8") as file:
            lines = file.readlines()
        lines = lines[17:] # skip header
        angstrom = '\u212B'.encode('utf-8')
        assert lines == [
            "RDF calculation:\n",
            "\n",
            f"    Number of bins: {rdf.n_bins}\n",
            f"    Bin width:      {rdf.delta_r} {angstrom}\n",
            f"    Minimum radius: {rdf.r_min} {angstrom}\n",
            f"    Maximum radius: {rdf.r_max} {angstrom}\n",
            "\n",
            f"    Number of frames: {rdf.n_frames}\n",
            f"    Number of atoms:  {rdf.n_atoms}\n",
            "\n",
            f"    Reference selection: {rdf.reference_selection}\n",
            f"    total number of atoms in reference selection: {len(rdf.reference_indices)}\n",
            f"    Target selection:    {rdf.target_selection}\n",
            f"    total number of atoms in target selection:    {len(rdf.target_indices)}\n",
            "\n",
            f"    Eliminate intra molecular contributions: {rdf.no_intra_molecular}\n",
            "\n",
            "\n",
            "\n", 
            "\n",
            "\n",
        ]


