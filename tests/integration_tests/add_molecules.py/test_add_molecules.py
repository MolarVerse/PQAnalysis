import pathlib
import pytest

from unittest.mock import patch
from filecmp import cmp as filecmp

pytestmark = pytest.mark.integration

from PQAnalysis.cli.add_molecules import main

__path__ = pathlib.Path(__file__).parent.absolute()



class TestAddMolecules:

    """
    Test the add_molecules function.
    """

    @pytest.mark.parametrize("example_dir", ["add_molecules"], indirect=False)
    def test_add_molecules_no_topology(self, test_integration_folder):
        """
        Test the add_molecules function when no topology is given.
        """
        # test command line tool with certain arguments

        with patch(
            'argparse._sys.argv',
            [
                'add_molecules',
                'mil68ga-md-01.rst',
                'perylene-md-05.xyz',
                '-o',
                'combined.rst'
            ]
        ):
            main()

        filecmp("combined.rst", str(__path__ / "ref_data" / "combined.rst"))

    @pytest.mark.parametrize("example_dir", ["add_molecules"], indirect=False)
    def test_add_molecules(self, test_integration_folder):
        """
        Test the add_molecules function when no topology is given.
        """
        # test command line tool with certain arguments

        with patch(
            'argparse._sys.argv',
            [
                'add_molecules',
                'mil68ga-md-01.rst',
                'perylene-md-05.xyz',
                '--top-file',
                'shake_mil.top',
                '--added-topology-file',
                'shake_perylene.top',
                '-o',
                'combined.rst',
                '--output-topology-file',
                'top.top'
            ]
        ):
            main()

        filecmp("combined.rst", str(__path__ / "ref_data" / "combined.rst"))
        filecmp("top.top", str(__path__ / "ref_data" / "top.top"))
