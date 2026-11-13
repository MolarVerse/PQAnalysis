"""
Integration test of the RDF between the centers of mass of molecules through
the input file and the rdf API on the 2000-molecule TIP3P water box
(examples/h2o_mm-mshake).

The result must equal the manual workflow that the keys replace: write the
center of mass trajectory with traj2com and run an atom RDF on it.
"""

import numpy as np
import pytest

from PQAnalysis.analysis.rdf import rdf
from PQAnalysis.io import rst2xyz
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis.tools import traj2com

pytestmark = pytest.mark.integration

DELTA_R = 0.1
R_MAX = 6.0


def write_input(name, **keys):
    lines = {
        "traj_files": "water.xyz",
        "reference_selection": "all",
        "target_selection": "all",
        "delta_r": DELTA_R,
        "r_max": R_MAX,
        "out_file": f"{name}.dat",
        "log_file": f"{name}.log",
    }
    lines.update(keys)
    lines = {key: value for key, value in lines.items() if value is not None}

    with open(f"{name}.in", "w", encoding="utf-8") as file:
        for key, value in lines.items():
            print(f"{key} = {value}", file=file)

    return f"{name}.in"


def read_table(name):
    return np.loadtxt(f"{name}.dat")


@pytest.mark.parametrize("example_dir", ["h2o_mm-mshake"], indirect=False)
class TestRDFCentersInputFile:

    @pytest.fixture(autouse=True)
    def water(self, test_integration_folder):
        rst2xyz("input.rst", output="water.xyz")

    def test_equals_the_manual_center_of_mass_workflow(self):
        # the manual workflow: a center of mass trajectory, then an atom RDF
        traj2com(
            ["water.xyz"],
            output="com.xyz",
            by_residue=True,
            restart_file="input.rst",
            moldescriptor_file="moldescriptor.dat",
        )
        rdf(
            write_input(
                "manual",
                traj_files="com.xyz",
                reference_selection="all",
                target_selection="all",
            )
        )

        rdf(
            write_input(
                "keys",
                restart_file="input.rst",
                moldescriptor_file="moldescriptor.dat",
                reference_com="True",
                target_com="True",
            )
        )

        manual, keys = read_table("manual"), read_table("keys")

        # exactly the same histogram, so exactly the same output
        assert np.array_equal(keys, manual)
        assert keys[:, 1].max() == pytest.approx(3.7, abs=0.4)

    def test_log_counts_molecules(self):
        rdf(
            write_input(
                "keys",
                restart_file="input.rst",
                moldescriptor_file="moldescriptor.dat",
                reference_com="True",
                target_com="True",
            )
        )

        log = open("keys.log", encoding="utf-8").read()

        assert (
            "total number of molecules (centers of mass) in reference "
            "selection: 2000" in log
        )

    def test_the_topology_files_are_inferred(self):
        # water.xyz -> water.rst, moldescriptor.dat in the same directory
        import shutil
        shutil.copy("input.rst", "water.rst")

        rdf(write_input("inferred", reference_com="True", target_com="True"))
        rdf(
            write_input(
                "explicit",
                restart_file="water.rst",
                moldescriptor_file="moldescriptor.dat",
                reference_com="True",
                target_com="True",
            )
        )

        assert np.array_equal(read_table("inferred"), read_table("explicit"))

    def test_one_side_with_and_without_the_intramolecular_pair(self):
        common = {
            "restart_file": "input.rst",
            "moldescriptor_file": "moldescriptor.dat",
            "reference_selection": "O",
            "target_com": "True",
            "reference_com": "False",
            "target_selection": "all",
        }
        rdf(write_input("kept", no_intra_molecular="False", **common))
        rdf(write_input("removed", no_intra_molecular="True", **common))
        rdf(write_input("default", **common))

        kept, removed, default = (
            read_table(name) for name in ("kept", "removed", "default")
        )

        # the oxygen is 0.065 A from the center of its own molecule: that
        # pair fills the first bin unless it is excluded (the default when
        # both topology files are given)
        assert kept[0, 1] > 100.0
        assert removed[0, 1] < 5.0
        assert np.array_equal(default, removed)

    def test_the_keys_need_the_topology_files(self):
        with pytest.raises(InputFileError, match="reference_com or target_com"):
            rdf(write_input("none", reference_com="True", target_com="True"))

    def test_without_the_keys_the_atom_rdf_is_unchanged(self):
        rdf(write_input("atoms"))
        rdf(
            write_input(
                "false", reference_com="False", target_com="False"
            )
        )

        assert np.array_equal(read_table("atoms"), read_table("false"))
