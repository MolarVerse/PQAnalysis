"""Tests for the ParmEd-based force-field import (parm2pq)."""

import sys

import numpy as np
import pytest

from PQAnalysis.io.force_field import (
    Parm2PQError,
    UnsupportedForceFieldTermError,
    convert_structure,
    load_structure,
    parm2pq,
)
from PQAnalysis.io.force_field.api import import_parmed

parmed = pytest.importorskip("parmed")

# pylint: disable=import-outside-toplevel,redefined-outer-name

_SCEE, _SCNB = 1.2, 2.0



def _atom_type(structure, name, rmin, epsilon, atomic_number, mass):
    from parmed.topologyobjects import AtomType
    types = getattr(structure, "_test_types", {})
    if name not in types:
        atom_type = AtomType(name, len(types) + 1, mass, atomic_number)
        atom_type.set_lj_params(epsilon, rmin)
        types[name] = atom_type
        structure._test_types = types  # pylint: disable=protected-access
    return types[name]



def _add_atom(structure, spec, residue, number, chain=""):
    from parmed.topologyobjects import Atom
    name, type_name, charge, rmin, epsilon, atomic_number, mass = spec
    atom = Atom(
        name=name,
        type=type_name,
        charge=charge,
        mass=mass,
        atomic_number=atomic_number,
        rmin=rmin,
        epsilon=epsilon,
    )
    atom.atom_type = _atom_type(
        structure, type_name, rmin, epsilon, atomic_number, mass
    )
    structure.add_atom(atom, residue, number, chain=chain)
    return atom



def _type(table, parm_type):
    table.append(parm_type)
    table.claim()
    return parm_type



def _add_chain(structure, n_atoms, residue, number, chain="", phase=30.0):
    """A linear chain C-C-...-C with one H at the first atom."""
    from parmed.topologyobjects import (
        Angle,
        AngleType,
        Bond,
        BondType,
        Dihedral,
        DihedralType,
    )
    atoms = [
        _add_atom(
            structure,
            ("C", "CT", 0.1 * (-1)**k, 1.9, 0.1, 6, 12.01),
            residue,
            number,
            chain,
        ) for k in range(n_atoms)
    ]
    hydrogen = _add_atom(
        structure,
        ("H", "HC", 0.05, 1.4, 0.02, 1, 1.008),
        residue,
        number,
        chain,
    )
    bond_type = _type(structure.bond_types, BondType(300.0, 1.5))
    h_type = _type(structure.bond_types, BondType(340.0, 1.09))
    angle_type = _type(structure.angle_types, AngleType(50.0, 109.5))
    dihedral_type = _type(
        structure.dihedral_types,
        DihedralType(1.4, 3, phase, scee=_SCEE, scnb=_SCNB),
    )
    for k in range(n_atoms - 1):
        structure.bonds.append(Bond(atoms[k], atoms[k + 1], bond_type))
    structure.bonds.append(Bond(atoms[0], hydrogen, h_type))
    for k in range(n_atoms - 2):
        structure.angles.append(
            Angle(atoms[k], atoms[k + 1], atoms[k + 2], angle_type)
        )
    for k in range(n_atoms - 3):
        structure.dihedrals.append(
            Dihedral(
                atoms[k],
                atoms[k + 1],
                atoms[k + 2],
                atoms[k + 3],
                type=dihedral_type,
            )
        )
    structure.dihedrals.append(
        Dihedral(hydrogen, atoms[0], atoms[1], atoms[2], type=dihedral_type)
    )
    return atoms, hydrogen



def _add_ion(structure, number):
    _add_atom(
        structure, ("NA", "Na+", 1.0, 1.37, 0.0028, 11, 22.99), "NA", number
    )



def _finish(structure):
    rng = np.random.default_rng(1)
    structure.coordinates = rng.uniform(0, 10, (len(structure.atoms), 3))
    return structure



@pytest.fixture
def chains():
    """Two identical 6-atom chains (each + one H) and a sodium ion."""
    structure = parmed.Structure()
    _add_chain(structure, 6, "BUT", 1)
    _add_ion(structure, 2)
    _add_chain(structure, 6, "BUT", 3)
    return _finish(structure)



def test_moltypes_and_atom_order(chains):
    setup = convert_structure(chains)

    assert [m.name for m in setup.moltypes] == ["BUT", "NA"]
    assert [m.n_atoms for m in setup.moltypes] == [7, 1]
    assert list(setup.atom_moltypes) == [1] * 7 + [2] + [1] * 7
    assert setup.moltypes[0].elements == ["C"] * 6 + ["H"]
    assert setup.moltypes[0].atom_types == [0] * 6 + [1]
    np.testing.assert_allclose(setup.positions, chains.coordinates)



def test_intra_pairs_exclusions_and_14(chains):
    pairs = convert_structure(chains).moltypes[0].intra_pairs

    # chain 1-2-3-4-5-6 and H7 bonded to 1
    assert pairs[0] == [-4, 5, 6]
    assert pairs[1] == [-5, 6]
    assert pairs[2] == [-6, -7]
    assert pairs[3] == [7]
    assert pairs[4] == [7]
    assert pairs[5] == [7]
    assert pairs[6] == []



def test_bonded_conventions(chains):
    setup = convert_structure(chains)

    assert setup.bond_types == [(1.5, 600.0), (1.09, 680.0)]
    assert setup.angle_types == [(109.5, 100.0)]
    assert setup.dihedral_types == [(1.4, 3, 330.0)]
    assert setup.scale_14_coulomb == pytest.approx(1 / _SCEE)
    assert setup.scale_14_vdw == pytest.approx(1 / _SCNB)
    assert len(setup.bonds) == 12
    assert setup.bonds[0] == (0, 1, 0, False)
    assert setup.dihedrals[0] == (0, 1, 2, 3, 0, False)
    assert not any(term[-1] for term in setup.bonds + setup.dihedrals)



def test_lennard_jones_pairs(chains):
    setup = convert_structure(chains)
    pairs = {(i, j): (c6, c12) for i, j, c6, c12 in setup.lj_pairs}

    assert set(pairs) == {(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)}
    rmin, epsilon = 1.9 + 1.4, np.sqrt(0.1 * 0.02)
    assert pairs[(0, 1)] == pytest.approx(
        (-2 * epsilon * rmin**6, epsilon * rmin**12)
    )



def test_shake_h_replaces_hydrogen_bonds(chains):
    setup = convert_structure(chains, shake_h=True)

    assert setup.shake == [(0, 6, 1.09), (8, 14, 1.09)]
    assert len(setup.bonds) == 10
    assert setup.bond_types == [(1.5, 600.0)]
    assert "shake = on;" in setup.input_lines(
        {
            "start": "s",
            "moldescriptor": "m",
            "topology": "t",
            "parameter": "p",
            "intra_nonbonded": "i",
        },
        12.0,
    )



def test_improper_goes_to_improper_section(chains):
    from parmed.topologyobjects import Dihedral, DihedralType
    atoms = chains.atoms
    improper_type = _type(
        chains.dihedral_types,
        DihedralType(10.5, 2, 180.0, scee=_SCEE, scnb=_SCNB)
    )
    chains.dihedrals.append(
        Dihedral(
            atoms[0],
            atoms[2],
            atoms[1],
            atoms[6],
            improper=True,
            ignore_end=True,
            type=improper_type,
        )
    )
    setup = convert_structure(chains)

    assert setup.improper_types == [(10.5, 2, 180.0)]
    assert setup.impropers == [(0, 2, 1, 6, 0)]
    assert "IMPROPERS" in setup.topology_lines()



def test_written_files(chains, tmp_path):
    setup = convert_structure(chains, rcoulomb=10.0, source="test")
    paths = setup.write(tmp_path, rcoulomb=10.0)

    assert sorted(p.name for p in paths.values()) == sorted(
        [
            "moldescriptor.dat",
            "topology.top",
            "parameter.prm",
            "intra_nonbonded.dat",
            "start.rst",
            "pq_input.in",
        ]
    )
    moldescriptor = paths["moldescriptor"].read_text()
    assert "BUT 7 0.05" in moldescriptor
    assert "NA 1 1" in moldescriptor

    parameter = paths["parameter"].read_text()
    assert "TYPES 0 0 0 0 0 0.833333333333 0.5" in parameter
    assert "\nNONCOULOMBICS LJ\n" in parameter
    assert "\nDIHEDRALS\n0 1.4 3 330\nEND" in parameter

    intra = paths["intra_nonbonded"].read_text().splitlines()
    assert intra[3:6] == ["1", "1 -4 5 6", "2 -5 6"]
    assert intra[-3:] == ["# NA", "2", "END"]

    start = paths["start"].read_text().splitlines()
    extent = np.ptp(chains.coordinates, axis=0).max()
    assert start[0] == f"Box {extent + 20:.12g} {extent + 20:.12g} " \
        f"{extent + 20:.12g} 90 90 90"
    assert start[1].split()[:3] == ["C", "1", "1"]
    assert start[8].split()[:3] == ["Na", "1", "2"]
    assert start[9].split()[:3] == ["C", "1", "1"]
    assert len(start[1].split()) == 12

    assert "rcoulomb = 10;" in paths["input"].read_text()

    with pytest.raises(FileExistsError, match="--overwrite"):
        setup.write(tmp_path, rcoulomb=10.0)
    setup.write(tmp_path, rcoulomb=10.0, overwrite=True)



def test_linker_between_chains():
    from parmed.topologyobjects import Angle, Bond, Dihedral
    structure = parmed.Structure()
    first, _ = _add_chain(structure, 4, "AAA", 1, chain="A")
    second, hydrogen = _add_chain(structure, 4, "BBB", 2, chain="B")
    bond_type = structure.bond_types[0]
    angle_type = structure.angle_types[0]
    dihedral_type = structure.dihedral_types[0]
    structure.bonds.append(Bond(first[3], second[0], bond_type))
    for atoms in (
        (first[2], first[3], second[0]),
        (first[3], second[0], second[1]),
        (first[3], second[0], hydrogen),
    ):
        structure.angles.append(Angle(*atoms, angle_type))
    for atoms in (
        (first[1], first[2], first[3], second[0]),
        (first[2], first[3], second[0], second[1]),
        (first[3], second[0], second[1], second[2]),
        (first[2], first[3], second[0], hydrogen),
    ):
        structure.dihedrals.append(Dihedral(*atoms, type=dihedral_type))
    setup = convert_structure(_finish(structure))

    # two molecules of one moltype
    assert len(setup.moltypes) == 1
    assert list(setup.atom_moltypes) == [1] * 10
    assert setup.bonds[-1] == (3, 5, 0, True)
    assert [a[-1] for a in setup.angles[-3:]] == [True] * 3
    assert [d[-1] for d in setup.dihedrals[-4:]] == [True] * 4
    assert sum(term[-1] for term in setup.bonds + setup.angles) == 4
    assert sum(term[-1] for term in setup.dihedrals) == 4
    assert any(line.endswith(" *") for line in setup.topology_lines())



def test_linker_without_carrier_term_raises():
    from parmed.topologyobjects import Bond
    structure = parmed.Structure()
    first, _ = _add_chain(structure, 4, "AAA", 1, chain="A")
    second, _ = _add_chain(structure, 4, "BBB", 2, chain="B")
    structure.bonds.append(Bond(first[3], second[0], structure.bond_types[0]))

    with pytest.raises(Parm2PQError, match="linker correction"):
        convert_structure(_finish(structure))



def test_unsupported_terms(chains):
    from parmed.topologyobjects import (
        BondType,
        Cmap,
        Improper,
        UreyBradley,
    )
    atoms = chains.atoms
    cases = [
        (
            "urey_bradleys",
            UreyBradley(atoms[0], atoms[2], BondType(1, 1)),
            "Urey-Bradley"
        ),
        ("impropers", Improper(*atoms[:4]), "harmonic"),
        ("cmaps", Cmap(*atoms[:5]), "CMAP"),
    ]
    for attribute, term, name in cases:
        getattr(chains, attribute).append(term)
        with pytest.raises(UnsupportedForceFieldTermError, match=name):
            convert_structure(chains)
        getattr(chains, attribute).clear()



def test_virtual_site_is_unsupported(chains):
    from parmed.topologyobjects import ExtraPoint
    chains.add_atom(ExtraPoint(name="EP", type="EP"), "BUT", 1)
    chains.coordinates = np.zeros((len(chains.atoms), 3))

    with pytest.raises(UnsupportedForceFieldTermError, match="virtual"):
        convert_structure(chains)



def test_special_14_lennard_jones_is_unsupported(chains):
    chains.atoms[0].epsilon_14 = 0.05

    with pytest.raises(UnsupportedForceFieldTermError, match="1-4 Lennard"):
        convert_structure(chains)



def test_non_uniform_14_scaling_is_unsupported(chains):
    chains.dihedrals[0].type.scee = 1.2
    from parmed.topologyobjects import DihedralType
    other = _type(
        chains.dihedral_types, DihedralType(1.0, 2, 0.0, scee=1.0, scnb=2.0)
    )
    chains.dihedrals[1].type = other

    with pytest.raises(UnsupportedForceFieldTermError, match="non-uniform"):
        convert_structure(chains)



def test_missing_coordinates_and_parameters(chains):
    chains.coordinates = None
    with pytest.raises(Parm2PQError, match="coordinates"):
        convert_structure(chains)

    _finish(chains)
    chains.bonds[0].type = None
    with pytest.raises(Parm2PQError, match="no parameters"):
        convert_structure(chains)



def test_amber_round_trip_and_cli(chains, tmp_path, monkeypatch):
    from parmed.amber import AmberParm
    prmtop, inpcrd = tmp_path / "x.prmtop", tmp_path / "x.rst7"
    chains.box = [30.0, 30.0, 30.0, 90.0, 90.0, 90.0]
    chains.velocities = np.full((len(chains.atoms), 3), 0.5)
    parm = AmberParm.from_structure(chains)
    parm.save(str(prmtop))
    parm.save(str(inpcrd))

    loaded = load_structure(prmtop, inpcrd)
    direct = convert_structure(chains)
    amber = convert_structure(loaded)

    np.testing.assert_allclose(amber.box, chains.box)
    np.testing.assert_allclose(amber.velocities, 0.5e12, rtol=1e-6)

    assert amber.scale_14_coulomb == pytest.approx(direct.scale_14_coulomb)
    np.testing.assert_allclose(
        np.array(amber.lj_pairs),
        np.array(direct.lj_pairs),
        rtol=1e-5,
    )
    assert amber.moltypes[0].intra_pairs == direct.moltypes[0].intra_pairs

    from PQAnalysis.cli.parm2pq import main
    out = tmp_path / "out"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "parm2pq",
            str(prmtop),
            "-c",
            str(inpcrd),
            "-d",
            str(out),
            "--shake-h",
            "--rcoulomb",
            "9",
        ],
    )
    main()
    assert (out / "topology.top").read_text().count("\nSHAKE\n") == 1
    assert "rcoulomb = 9;" in (out / "pq_input.in").read_text()

    setup = parm2pq(prmtop, inpcrd, output_dir=out, overwrite=True)
    assert setup.source == "x.prmtop + x.rst7"



def test_charmm_round_trip(tmp_path, chains):
    from parmed.charmm import CharmmParameterSet
    psf, crd = tmp_path / "x.psf", tmp_path / "x.crd"
    prm, rtf = tmp_path / "x.prm", tmp_path / "x.rtf"
    CharmmParameterSet.from_structure(chains).write(par=str(prm), top=str(rtf))
    chains.save(str(psf))
    chains.save(str(crd))

    with pytest.raises(Parm2PQError, match="parameter files"):
        load_structure(psf)

    direct = convert_structure(chains)
    charmm = convert_structure(load_structure(psf, crd, [rtf, prm]))

    # CHARMM stores the 1-4 Lennard-Jones depth instead of SCNB
    assert charmm.scale_14_coulomb == pytest.approx(direct.scale_14_coulomb)
    assert charmm.scale_14_vdw == pytest.approx(direct.scale_14_vdw)
    assert charmm.moltypes[0].intra_pairs == direct.moltypes[0].intra_pairs
    np.testing.assert_allclose(
        np.array(charmm.lj_pairs),
        np.array(direct.lj_pairs),
        rtol=1e-6,
    )
    assert charmm.dihedral_types == direct.dihedral_types
    np.testing.assert_allclose(charmm.positions, direct.positions, atol=1e-4)



def test_missing_parmed_message(monkeypatch):
    monkeypatch.setitem(sys.modules, "parmed", None)
    with pytest.raises(ImportError, match=r"PQAnalysis\[ff\]"):
        import_parmed()



def _adjusts(structure, epsilon_scale, chgscale=0.5, rmin_shift=0.0):
    from parmed.topologyobjects import (
        NonbondedException,
        NonbondedExceptionType,
    )
    for dihedral in structure.dihedrals:
        atom_1, atom_4 = dihedral.atom1, dihedral.atom4
        rmin = atom_1.rmin + atom_4.rmin + rmin_shift
        epsilon = np.sqrt(atom_1.epsilon * atom_4.epsilon) * epsilon_scale
        exception_type = _type(
            structure.adjust_types,
            NonbondedExceptionType(rmin, epsilon, chgscale),
        )
        structure.adjusts.append(
            NonbondedException(atom_1, atom_4, type=exception_type)
        )



def test_adjusts_set_14_scaling(chains):
    _adjusts(chains, 0.5)
    setup = convert_structure(chains)

    assert setup.scale_14_coulomb == pytest.approx(0.5)
    assert setup.scale_14_vdw == pytest.approx(0.5)
    assert setup.moltypes[0].intra_pairs[0] == [-4, 5, 6]



def test_adjusts_with_own_lennard_jones_are_unsupported(chains):
    _adjusts(chains, 0.5, rmin_shift=0.3)
    with pytest.raises(UnsupportedForceFieldTermError, match="pair-specific"):
        convert_structure(chains)



def test_non_uniform_adjusts_are_unsupported(chains):
    _adjusts(chains, 0.5)
    chains.adjusts[0].type.chgscale = 0.8
    with pytest.raises(UnsupportedForceFieldTermError, match="non-uniform"):
        convert_structure(chains)



def test_nbfix_and_geometric_rule(chains):
    chains.atoms[0].atom_type.nbfix = {"HC": (3.0, 0.2, 3.0, 0.2)}
    pairs = {
        (i, j): (c6, c12)
        for i, j, c6, c12 in convert_structure(chains).lj_pairs
    }
    assert pairs[(0, 1)] == pytest.approx((-2 * 0.2 * 3.0**6, 0.2 * 3.0**12))

    chains.atoms[0].atom_type.nbfix = {}
    chains.combining_rule = "geometric"
    pairs = {
        (i, j): (c6, c12)
        for i, j, c6, c12 in convert_structure(chains).lj_pairs
    }
    atom_c, atom_h = chains.atoms[0], chains.atoms[6]
    sigma = np.sqrt(atom_c.sigma * atom_h.sigma)
    epsilon = np.sqrt(0.1 * 0.02)
    assert pairs[(0, 1)] == pytest.approx(
        (-4 * epsilon * sigma**6, 4 * epsilon * sigma**12)
    )



def test_dihedral_without_parameters_and_hbond_terms(chains):
    chains.parm_data = {"HBOND_ACOEF": [1.0]}
    with pytest.raises(UnsupportedForceFieldTermError, match="10-12"):
        convert_structure(chains)

    chains.parm_data = {}
    chains.dihedrals[0].type = None
    with pytest.raises(Parm2PQError, match="no parameters"):
        convert_structure(chains)



def test_periodic_box_and_velocities(chains):
    chains.box = [30.0, 31.0, 32.0, 90.0, 90.0, 90.0]
    chains.velocities = np.ones((len(chains.atoms), 3))
    setup = convert_structure(chains)

    np.testing.assert_allclose(setup.box, chains.box)
    np.testing.assert_allclose(setup.velocities, 1e12)
    assert setup.restart_lines()[0] == "Box 30 31 32 90 90 90"
    assert setup.restart_lines()[1].split()[6] == "1.00000000e+12"



def test_special_14_radius_is_unsupported(chains):
    chains.atoms[7].rmin_14 = 2.0

    with pytest.raises(UnsupportedForceFieldTermError, match="1-4 Lennard"):
        convert_structure(chains)



def test_zero_14_scaling_constant_raises(chains):
    chains.dihedrals[0].type.scee = 0.0

    with pytest.raises(Parm2PQError, match="scaling constant of 0"):
        convert_structure(chains)



def test_adjust_without_parameters_raises(chains):
    _adjusts(chains, 0.5)
    chains.adjusts[0].type = None

    with pytest.raises(Parm2PQError, match="1-4 pair has no parameters"):
        convert_structure(chains)



def test_adjusts_skip_pairs_without_lennard_jones(chains):
    for atom in chains.atoms:
        if atom.atomic_number == 1:
            atom.epsilon = atom.epsilon_14 = 0.0
    _adjusts(chains, 0.5)
    setup = convert_structure(chains)

    assert setup.scale_14_vdw == pytest.approx(0.5)



def test_angle_without_parameters_raises(chains):
    chains.angles[0].type = None

    with pytest.raises(Parm2PQError, match="angle 1-2-3 has no parameters"):
        convert_structure(chains)



def test_amber_pair_without_lennard_jones(chains):
    from parmed.amber import AmberParm
    parm = AmberParm.from_structure(chains)
    n_types = parm.ptr("ntypes")
    nb_idx = parm.atoms[0].nb_idx
    parm.parm_data["NONBONDED_PARM_INDEX"][n_types * (nb_idx - 1) + nb_idx -
                                           1] = 0
    setup = convert_structure(parm)

    assert setup.lj_pairs[0] == (0, 0, 0.0, 0.0)



def test_moltype_names_are_unique_and_not_numeric():
    structure = parmed.Structure()
    _add_chain(structure, 4, "AAA", 1)
    _add_chain(structure, 5, "AAA", 2)
    _add_atom(structure, ("CL", "Cl-", -1.0, 2.47, 0.1, 17, 35.45), "7", 3)
    setup = convert_structure(_finish(structure))

    assert [m.name for m in setup.moltypes] == ["AAA_1", "AAA_2", "MOL7"]



def test_restart_counter_restarts_for_adjacent_molecules():
    structure = parmed.Structure()
    _add_chain(structure, 4, "AAA", 1)
    _add_chain(structure, 4, "AAA", 2)
    lines = convert_structure(_finish(structure)).restart_lines()

    assert [line.split()[1] for line in lines[1:]] == [
        str(k) for k in (1, 2, 3, 4, 5) * 2
    ]
