"""Crystal preparation and restart round-trip tests."""

import numpy as np
import pytest

from PQAnalysis.io.crystal import prepare_crystal
from PQAnalysis.io.restart_file.api import read_restart_file, write_restart_file

CIF = """data_example
_cell_length_a 8
_cell_length_b 9
_cell_length_c 10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma {gamma}
_space_group_name_H-M_alt 'P -1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
C1 C 0.2 0.3 0.4 {occupancy}
"""



@pytest.mark.parametrize("gamma", [90, 73])
def test_supercell_restart_round_trip(tmp_path, gamma):
    source = tmp_path / "crystal.cif"
    source.write_text(CIF.format(gamma=gamma, occupancy=1), encoding="ascii")
    system = prepare_crystal(str(source), (2, 1, 1))

    assert system.n_atoms == 4
    restart = str(tmp_path / "crystal.rst")
    write_restart_file(system, restart)
    restored = read_restart_file(restart)
    assert restored.n_atoms == 4
    np.testing.assert_allclose(restored.pos, system.pos, atol=1e-8)
    assert restored.cell.gamma == pytest.approx(gamma)
    assert restored.cell.x == pytest.approx(16)



def test_requires_ordered_sites_and_checks_atom_budget(tmp_path):
    source = tmp_path / "crystal.cif"
    source.write_text(CIF.format(gamma=90, occupancy=0.5), encoding="ascii")
    with pytest.raises(ValueError, match="partial occupancies"):
        prepare_crystal(str(source))
    source.write_text(CIF.format(gamma=90, occupancy=1), encoding="ascii")
    with pytest.raises(ValueError, match="atom limit"):
        prepare_crystal(str(source), (2, 1, 1), max_atoms=3)
    with pytest.raises(ValueError, match="positive integer repeats"):
        prepare_crystal(str(source), (0, 1, 1))



def test_requires_explicit_block_when_ambiguous(tmp_path):
    first = CIF.format(gamma=90, occupancy=1)
    second = CIF.replace("data_example",
                         "data_second").format(gamma=73, occupancy=1)
    source = tmp_path / "crystal.cif"
    source.write_text(first + second, encoding="ascii")
    with pytest.raises(ValueError, match="exactly one"):
        prepare_crystal(str(source))
    assert prepare_crystal(str(source), block="second").cell.gamma == 73



def test_cli_writes_restart(tmp_path, monkeypatch):
    from PQAnalysis.cli.main import main

    source = tmp_path / "crystal.cif"
    source.write_text(CIF.format(gamma=73, occupancy=1), encoding="ascii")
    output = tmp_path / "supercell.rst"
    monkeypatch.setattr(
        "sys.argv",
        [
            "pqanalysis",
            "crystal",
            str(source),
            str(output),
            "--repeat",
            "2",
            "1",
            "1",
        ]
    )
    main()
    assert read_restart_file(str(output)).n_atoms == 4



def test_rejects_equivalent_source_sites(tmp_path):
    source = tmp_path / "crystal.cif"
    source.write_text(
        CIF.format(gamma=90, occupancy=1) + "C2 C 0.8 0.7 0.6 1\n",
        encoding="ascii",
    )
    with pytest.raises(ValueError, match="equivalent"):
        prepare_crystal(str(source))



def test_rejects_invalid_atom_limit():
    with pytest.raises(ValueError, match="atom limit must be positive"):
        prepare_crystal("unused.cif", max_atoms=0)



@pytest.mark.parametrize(
    ("original", "replacement", "message"),
    [
        ("_cell_length_a 8", "_cell_length_a nan", "finite, nonsingular"),
        ("0.2 0.3 0.4", "1e309 0.3 0.4", "finite atomic positions"),
    ],
)
def test_rejects_invalid_geometry(tmp_path, original, replacement, message):
    source = tmp_path / "crystal.cif"
    source.write_text(
        CIF.format(gamma=90, occupancy=1).replace(original, replacement),
        encoding="ascii",
    )
    with pytest.raises(ValueError, match=message):
        prepare_crystal(str(source))



def test_rejects_near_coincident_periodic_sites(tmp_path):
    source = tmp_path / "crystal.cif"
    source.write_text(
        CIF.format(gamma=90, occupancy=1).replace("P -1", "P 1")
        + "C2 C 1.2011 0.3 0.4 1\n",
        encoding="ascii",
    )
    with pytest.raises(ValueError, match="coincident periodic sites"):
        prepare_crystal(str(source))



def test_cli_reports_invalid_cif(tmp_path, monkeypatch):
    from PQAnalysis.cli.main import main

    source = tmp_path / "crystal.cif"
    source.write_text(CIF.format(gamma=90, occupancy=0.5), encoding="ascii")
    output = tmp_path / "supercell.rst"
    monkeypatch.setattr(
        "sys.argv", ["pqanalysis", "crystal", str(source), str(output)]
    )
    with pytest.raises(SystemExit, match="partial occupancies"):
        main()
    assert not output.exists()
