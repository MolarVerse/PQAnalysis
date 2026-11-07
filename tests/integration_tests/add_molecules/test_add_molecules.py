import pathlib
import pytest

import numpy as np

from unittest.mock import patch
from filecmp import cmp as filecmp

pytestmark = pytest.mark.integration

from PQAnalysis.cli.add_molecules import main

__path__ = pathlib.Path(__file__).parent.absolute()

EXAMPLES = __path__.parents[2] / "examples" / "add_molecules"
CUTOFF = 1.0  # default of the --cut option in Angstrom


def read_restart(path):
    """
    Splits an orthorhombic restart file into its box lengths and atom rows.

    Each row is a list of tokens: name, two integer columns and the floats
    (x y z, optionally followed by velocities and forces).
    """

    rows = []
    box = None

    for line in pathlib.Path(path).read_text().splitlines():
        tokens = line.split()

        if not tokens or tokens[0] == "Step":
            continue
        if tokens[0] == "Box":
            box = np.array([float(value) for value in tokens[1:4]])
            continue

        rows.append(tokens)

    return box, rows


def positions(rows):
    return np.array([[float(value) for value in row[3:6]] for row in rows])


def distance_matrix(pos, box):
    diff = pos[:, None, :] - pos[None, :, :]
    diff -= box * np.round(diff / box)

    return np.linalg.norm(diff, axis=-1)


def assert_added_molecule(combined, reference, molecule):
    """
    Asserts that the combined restart file is the reference up to the
    placement of the added molecule.

    The framework and the box must be identical to the reference. The
    placement depends on the random number stream (seeded by the tests),
    which can change between numpy/scipy versions, so for the added atoms
    only what the algorithm guarantees is checked: names and columns, the
    shape of the molecule (rigid placement) and the distance cutoff to the
    framework.
    """

    box, rows = read_restart(combined)
    ref_box, ref_rows = read_restart(reference)
    n_added = len(molecule)

    assert np.array_equal(box, ref_box)
    assert len(rows) == len(ref_rows)

    framework, added = rows[:-n_added], rows[-n_added:]
    ref_framework, ref_added = ref_rows[:-n_added], ref_rows[-n_added:]

    assert framework == ref_framework
    assert [row[:3] for row in added] == [row[:3] for row in ref_added]

    pos = positions(added)
    assert np.allclose(
        distance_matrix(pos, box), distance_matrix(molecule, box), atol=1e-5
    )

    diff = pos[:, None, :] - positions(framework)[None, :, :]
    diff -= box * np.round(diff / box)
    assert np.linalg.norm(diff, axis=-1).min() >= CUTOFF - 1e-6


def read_xyz_positions(path):
    lines = pathlib.Path(path).read_text().splitlines()[2:]

    return np.array([[float(value) for value in line.split()[1:4]]
                     for line in lines if line.strip()])



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

        assert_added_molecule(
            "combined.rst",
            __path__ / "ref_data" / "combined.rst",
            read_xyz_positions("perylene-md-05.xyz"),
        )

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

        assert_added_molecule(
            "combined.rst",
            __path__ / "ref_data" / "combined.rst",
            read_xyz_positions("perylene-md-05.xyz"),
        )
        assert filecmp("top.top", str(__path__ / "ref_data" / "top.top"))



class TestAssertAddedMolecule:

    """
    The comparison must accept the reference and reject every change that
    matters, otherwise the integration tests above would prove nothing.
    """

    REFERENCE = __path__ / "ref_data" / "combined.rst"

    @pytest.fixture
    def molecule(self):
        return read_xyz_positions(EXAMPLES / "perylene-md-05.xyz")

    def modified(self, tmp_path, edit):
        """Writes a copy of the reference whose atom lines are edited."""

        lines = self.REFERENCE.read_text().splitlines()
        atom_lines = [
            i for i, line in enumerate(lines)
            if line.split() and line.split()[0] not in ("Box", "Step")
        ]
        edit(lines, atom_lines)

        path = tmp_path / "modified.rst"
        path.write_text("\n".join(lines) + "\n")

        return path

    @staticmethod
    def shifted(line, delta):
        tokens = line.split()
        tokens[3:6] = [str(float(value) + d) for value, d in zip(tokens[3:6], delta)]

        return " ".join(tokens)

    def test_accepts_the_reference(self, molecule):
        assert_added_molecule(self.REFERENCE, self.REFERENCE, molecule)

    def test_rejects_a_changed_framework_atom(self, tmp_path, molecule):
        def edit(lines, atoms):
            lines[atoms[0]] = self.shifted(lines[atoms[0]], (0.01, 0, 0))

        with pytest.raises(AssertionError):
            assert_added_molecule(
                self.modified(tmp_path, edit), self.REFERENCE, molecule
            )

    def test_rejects_a_distorted_added_molecule(self, tmp_path, molecule):
        def edit(lines, atoms):
            lines[atoms[-1]] = self.shifted(lines[atoms[-1]], (0.5, 0, 0))

        with pytest.raises(AssertionError):
            assert_added_molecule(
                self.modified(tmp_path, edit), self.REFERENCE, molecule
            )

    def test_rejects_an_added_molecule_overlapping_the_framework(
        self, tmp_path, molecule
    ):
        _, rows = read_restart(self.REFERENCE)
        n_added = len(molecule)
        target = positions(rows)[0] + 0.1
        delta = target - positions(rows[-n_added:])[0]

        def edit(lines, atoms):
            for index in atoms[-n_added:]:
                lines[index] = self.shifted(lines[index], delta)

        with pytest.raises(AssertionError):
            assert_added_molecule(
                self.modified(tmp_path, edit), self.REFERENCE, molecule
            )

    def test_rejects_a_missing_atom(self, tmp_path, molecule):
        def edit(lines, atoms):
            del lines[atoms[-1]]

        with pytest.raises(AssertionError):
            assert_added_molecule(
                self.modified(tmp_path, edit), self.REFERENCE, molecule
            )
