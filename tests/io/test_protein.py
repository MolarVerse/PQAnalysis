"""Regression tests for protein preparation."""

from math import isfinite
from pathlib import Path
import sys


import pytest

from PQAnalysis.io.protein import prepare_protein



@pytest.mark.parametrize("padding", [0.0, float("nan"), float("inf")])
def test_invalid_padding(tmp_path, padding):
    """Reject invalid padding before reading input."""
    with pytest.raises(ValueError, match="padding"):
        prepare_protein(tmp_path / "missing.pdb", tmp_path, padding_nm=padding)



@pytest.mark.parametrize("filename", ["system.prmtop", "system.inpcrd"])
def test_existing_output(tmp_path, filename):
    """Never overwrite either AMBER file."""
    existing = tmp_path / filename
    existing.write_text("keep")
    with pytest.raises(FileExistsError):
        prepare_protein(tmp_path / "missing.pdb", tmp_path)
    assert existing.read_text() == "keep"
    assert list(tmp_path.iterdir()) == [existing]



@pytest.mark.parametrize("residue", ["UNK", "ALA"])
def test_unparameterizable_residue(tmp_path, residue):
    """Reject an unknown residue or alanine missing heavy atoms."""
    pytest.importorskip("openmm")
    pytest.importorskip("parmed")
    pdb = tmp_path / "input.pdb"
    pdb.write_text(
        f"ATOM      1  CA  {residue} A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    with pytest.raises(ValueError, match="template"):
        prepare_protein(pdb, tmp_path)
    assert list(tmp_path.iterdir()) == [pdb]


def test_prepare_protein_round_trip(tmp_path):
    """Solvate OpenMM's bundled protein and reopen the AMBER pair."""
    app = pytest.importorskip("openmm.app")
    parmed = pytest.importorskip("parmed")
    fixture = Path(app.__file__).with_name("data") / "test.pdb"
    pdb = tmp_path / "protein.pdb"
    pdb.write_text("".join(
        line for line in fixture.read_text().splitlines(keepends=True)
        if line.startswith("ATOM  ") and line[17:20].strip() != "HOH"
    ) + "END\n")

    result = prepare_protein(pdb, tmp_path)
    system = parmed.load_file(str(result.topology), xyz=str(result.coordinates))
    assert result.water_molecules > 0
    assert result.total_atoms == len(system.atoms)
    assert result.box_angstrom == pytest.approx(system.box[:3])
    assert isfinite(result.energy_kcal_mol)


def test_missing_openmm(tmp_path, monkeypatch):
    """Explain how to install the optional preparation dependencies."""
    monkeypatch.setitem(sys.modules, "openmm", None)
    with pytest.raises(ImportError, match=r"PQAnalysis\[protein\]"):
        prepare_protein(tmp_path / "protein.pdb", tmp_path)
