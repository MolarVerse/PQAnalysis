"""Regression tests for protein preparation."""

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
