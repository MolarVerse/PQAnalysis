"""
Conversion of a ParmEd structure into a :py:class:`PQForceFieldSetup`.

ParmEd stores every force field in AMBER units and conventions; the
conversion maps those onto the PQ energy expressions:

- bonds and angles: PQ uses 1/2 k, AMBER k, so k is doubled;
- proper and periodic improper torsions: PQ evaluates
  k (1 + cos(n phi + delta)) with the IUPAC sign of phi, AMBER
  k (1 + cos(n phi - delta)), so the phase changes sign;
- Lennard-Jones: one c6/c12 pair per pair of global van der Waals types;
- exclusions and scaled 1-4 pairs: the PQ intra non-bonded file inside a
  molecule, linker marks for pairs between molecules;
- rigid three-site water (O-H, O-H and H-H bond, as tleap writes TIP3P):
  with SHAKE, one M-SHAKE reference geometry per moltype instead of three
  coupled SHAKE constraints.
"""

import math

from collections import defaultdict, deque
from numbers import Real

import numpy as np

from beartype.typing import Any

from .exceptions import Parm2PQError, UnsupportedForceFieldTermError
from .pq_setup import PQForceFieldSetup, PQMolType

_REL_TOL = 1e-6

#: ParmEd structure lists that hold terms PQ cannot evaluate
_UNSUPPORTED_LISTS = (
    ("cmaps", "CHARMM CMAP correction maps"),
    ("urey_bradleys", "Urey-Bradley terms"),
    ("impropers", "harmonic (CHARMM) impropers"),
    ("rb_torsions", "Ryckaert-Bellemans torsions"),
    ("trigonal_angles", "trigonal angles (AMOEBA)"),
    ("out_of_plane_bends", "out-of-plane bends (AMOEBA)"),
    ("pi_torsions", "pi torsions (AMOEBA)"),
    ("stretch_bends", "stretch-bend terms (AMOEBA)"),
    ("torsion_torsions", "torsion-torsion terms (AMOEBA)"),
    ("multipole_frames", "multipoles (AMOEBA)"),
)



def _close(value_1: Real, value_2: Real) -> bool:
    return math.isclose(value_1, value_2, rel_tol=_REL_TOL, abs_tol=1e-10)



def check_supported(structure: Any) -> None:
    """
    Raises if the structure contains a term PQ cannot evaluate.

    Parameters
    ----------
    structure : parmed.Structure
        The structure to check.

    Raises
    ------
    UnsupportedForceFieldTermError
        For CMAP, Urey-Bradley, harmonic impropers, AMOEBA terms, virtual
        sites or extra points, special 1-4 Lennard-Jones parameters and
        AMBER 10-12 hydrogen-bond terms.
    """
    for attribute, name in _UNSUPPORTED_LISTS:
        terms = getattr(structure, attribute, None)
        if terms:
            raise UnsupportedForceFieldTermError(
                f"The force field contains {len(terms)} {name}; PQ cannot "
                "evaluate this term."
            )

    from parmed.topologyobjects import ExtraPoint  # pylint: disable=import-outside-toplevel

    virtual = [
        atom for atom in structure.atoms
        if isinstance(atom, ExtraPoint) or atom.atomic_number == 0
    ]
    if virtual:
        raise UnsupportedForceFieldTermError(
            f"The force field contains {len(virtual)} virtual sites or extra "
            f"points (first: atom {virtual[0].idx + 1} '{virtual[0].name}'); "
            "PQ cannot evaluate them."
        )

    lj_14_factor(structure)

    parm_data = getattr(structure, "parm_data", {})
    if any(abs(value) > 0 for value in parm_data.get("HBOND_ACOEF", [])):
        raise UnsupportedForceFieldTermError(
            "The force field contains AMBER 10-12 hydrogen-bond terms; PQ "
            "cannot evaluate them."
        )



def lj_14_factor(structure: Any) -> Real:
    """
    Returns the factor between the 1-4 and the normal Lennard-Jones
    well depth of the atoms (CHARMM files store the scaled 1-4 values).

    Parameters
    ----------
    structure : parmed.Structure
        The structure.

    Returns
    -------
    float
        The common factor epsilon_14 / epsilon, 1.0 without 1-4 values.

    Raises
    ------
    UnsupportedForceFieldTermError
        If the 1-4 radius differs or the factor is not the same for all
        atoms.
    """
    factors = {}
    for atom in structure.atoms:
        same_radius = _close(atom.rmin_14, atom.rmin)
        if atom.epsilon != 0.0 and same_radius:
            factor = atom.epsilon_14 / atom.epsilon
            factors.setdefault(round(factor, 6), factor)
        elif not same_radius or atom.epsilon_14 != 0.0:
            factors[None] = None
        if None in factors or len(factors) > 1:
            raise UnsupportedForceFieldTermError(
                "The force field has special 1-4 Lennard-Jones parameters "
                f"(atom {atom.idx + 1} '{atom.name}', type '{atom.type}'); "
                "PQ only scales the normal Lennard-Jones pair."
            )
    return next(iter(factors.values()), 1.0)



def _dihedral_terms(dihedral: Any) -> list:
    if dihedral.type is None:
        raise Parm2PQError(
            "The dihedral "
            f"{'-'.join(str(a.idx + 1) for a in _atoms_of(dihedral))} has no "
            "parameters."
        )
    try:
        return list(dihedral.type)
    except TypeError:
        return [dihedral.type]



def _atoms_of(term: Any) -> list:
    names = ("atom1", "atom2", "atom3", "atom4")
    return [getattr(term, name) for name in names if hasattr(term, name)]



def _combined_lj(atom_1: Any, atom_2: Any, rule: str) -> tuple[Real, Real]:
    """Returns (Rmin, epsilon) of an atom pair without NBFIX."""
    epsilon = math.sqrt(abs(atom_1.epsilon * atom_2.epsilon))
    if rule == "geometric":
        sigma = math.sqrt(atom_1.sigma * atom_2.sigma)
        return sigma * 2**(1 / 6), epsilon
    return atom_1.rmin + atom_2.rmin, epsilon



class _LennardJones:

    """Global van der Waals types and their c6/c12 pair parameters."""

    def __init__(self, structure: Any):
        self.structure = structure
        parm_data = getattr(structure, "parm_data", {})
        self.amber_tables = "LENNARD_JONES_ACOEF" in parm_data
        self.rule = getattr(structure, "combining_rule", "lorentz")
        self.keys = []
        self.representatives = []
        self.atom_vdw_types = []

        index_of = {}
        for atom in structure.atoms:
            key = self._key(atom)
            if key not in index_of:
                index_of[key] = len(self.keys)
                self.keys.append(key)
                self.representatives.append(atom)
            self.atom_vdw_types.append(index_of[key])

    def _key(self, atom: Any) -> tuple:
        if self.amber_tables:
            return (atom.nb_idx, )
        return (
            str(atom.atom_type), round(atom.rmin, 8), round(atom.epsilon, 8)
        )

    def pair(self, type_1: int, type_2: int) -> tuple[Real, Real]:
        """Returns (c6, c12) of two global van der Waals types."""
        atom_1 = self.representatives[type_1]
        atom_2 = self.representatives[type_2]
        if self.amber_tables:
            data = self.structure.parm_data
            n_types = self.structure.ptr("ntypes")
            index = data["NONBONDED_PARM_INDEX"][n_types *
                                                 (atom_1.nb_idx - 1) +
                                                 atom_2.nb_idx - 1]
            if index <= 0:
                return 0.0, 0.0
            return (
                -data["LENNARD_JONES_BCOEF"][index - 1],
                data["LENNARD_JONES_ACOEF"][index - 1],
            )

        nbfix = getattr(atom_1.atom_type, "nbfix", {}) or {}
        if str(atom_2.atom_type) in nbfix:
            rmin, epsilon = nbfix[str(atom_2.atom_type)][:2]
        else:
            rmin, epsilon = _combined_lj(atom_1, atom_2, self.rule)
        return -2 * epsilon * rmin**6, epsilon * rmin**12

    def all_pairs(self) -> list[tuple[int, int, Real, Real]]:
        """Returns all pairs i <= j with their c6 and c12."""
        n_types = len(self.keys)
        return [
            (i, j, *self.pair(i, j))
            for i in range(n_types)
            for j in range(i, n_types)
        ]



def _scale_from_dihedrals(
    structure: Any,
    lj: _LennardJones,
) -> tuple[Real, Real]:
    """Uniform 1-4 scaling; pairs without charge or LJ do not count."""
    scee_values, scnb_values = {}, {}
    for dihedral in structure.dihedrals:
        if dihedral.ignore_end:
            continue
        scee = scnb = 0.0
        for term in _dihedral_terms(dihedral):
            if scee == 0.0 or scnb == 0.0:
                scee, scnb = term.scee, term.scnb
        if scee == 0.0 or scnb == 0.0:
            raise Parm2PQError(
                "A dihedral with 1-4 interactions has a 1-4 scaling "
                "constant of 0."
            )
        atom_1, atom_4 = dihedral.atom1, dihedral.atom4
        if atom_1.charge * atom_4.charge != 0.0:
            scee_values.setdefault(round(scee, 6), scee)
        c6, c12 = lj.pair(
            lj.atom_vdw_types[atom_1.idx],
            lj.atom_vdw_types[atom_4.idx],
        )
        if c6 != 0.0 or c12 != 0.0:
            scnb_values.setdefault(round(scnb, 6), scnb)

    if len(scee_values) > 1 or len(scnb_values) > 1:
        raise UnsupportedForceFieldTermError(
            "The force field uses non-uniform 1-4 scaling (SCEE "
            f"{sorted(scee_values)}, SCNB {sorted(scnb_values)}); PQ has one "
            "global 1-4 scaling factor for Coulomb and one for Lennard-Jones."
        )
    scee = next(iter(scee_values.values()), 1.0)
    scnb = next(iter(scnb_values.values()), 1.0)
    return 1.0 / scee, 1.0 / scnb



def _scale_from_adjusts(structure: Any, rule: str) -> tuple[Real, Real]:
    coulomb, vdw = {}, {}
    for adjust in structure.adjusts:
        if adjust.type is None:
            raise Parm2PQError("A 1-4 pair has no parameters.")
        chgscale = adjust.type.chgscale
        coulomb.setdefault(round(chgscale, 6), chgscale)
        rmin, epsilon = _combined_lj(adjust.atom1, adjust.atom2, rule)
        if epsilon == 0.0 and adjust.type.epsilon == 0.0:
            continue
        if epsilon == 0.0 or not _close(adjust.type.rmin, rmin):
            raise UnsupportedForceFieldTermError(
                "The force field has pair-specific 1-4 Lennard-Jones "
                f"parameters (atoms {adjust.atom1.idx + 1} and "
                f"{adjust.atom2.idx + 1}); PQ only scales the normal "
                "Lennard-Jones pair."
            )
        ratio = adjust.type.epsilon / epsilon
        vdw.setdefault(round(ratio, 6), ratio)

    if len(coulomb) > 1 or len(vdw) > 1:
        raise UnsupportedForceFieldTermError(
            "The force field uses non-uniform 1-4 scaling (Coulomb "
            f"{sorted(coulomb)}, Lennard-Jones {sorted(vdw)}); PQ has one "
            "global 1-4 scaling factor for each."
        )
    return next(iter(coulomb.values()), 1.0), next(iter(vdw.values()), 1.0)



def _pair(index_1: int, index_2: int) -> tuple[int, int]:
    return (index_1, index_2) if index_1 < index_2 else (index_2, index_1)



def _nonbonded_pairs(structure: Any) -> tuple[set, set]:
    """
    Returns the fully excluded and the scaled 1-4 pairs as ParmEd and
    OpenMM build them: pairs up to nrexcl - 1 bonds apart and the ends of
    proper dihedrals are excluded, the ends of dihedrals that keep their
    1-4 (or the explicit 1-4 adjust list) are scaled, explicit exclusions
    win over everything.
    """
    n_excl = (getattr(structure, "nrexcl", 3) or 3) - 1
    excluded = set()
    for atom in structure.atoms:
        seen = {atom.idx: 0}
        queue = deque([atom])
        while queue:
            current = queue.popleft()
            if seen[current.idx] == n_excl:
                continue
            for partner in current.bond_partners:
                if partner.idx not in seen:
                    seen[partner.idx] = seen[current.idx] + 1
                    queue.append(partner)
        excluded.update(
            _pair(atom.idx, other) for other in seen if other != atom.idx
        )
        excluded.update(
            _pair(atom.idx, other.idx) for other in atom.exclusion_partners
        )

    if structure.adjusts:
        pairs_14 = {
            _pair(adjust.atom1.idx, adjust.atom2.idx)
            for adjust in structure.adjusts
        }
    else:
        excluded.update(
            _pair(dihedral.atom1.idx, dihedral.atom4.idx)
            for dihedral in structure.dihedrals
            if not dihedral.improper
        )
        pairs_14 = {
            _pair(dihedral.atom1.idx, dihedral.atom4.idx)
            for dihedral in structure.dihedrals
            if not dihedral.ignore_end
        }
    explicit = {
        _pair(atom.idx, other.idx)
        for atom in structure.atoms
        for other in atom.exclusion_partners
    }
    pairs_14 -= explicit
    excluded -= pairs_14
    return excluded, pairs_14



def _molecules(structure: Any) -> list[list[int]]:
    """Connected components, split by chain, sorted by first atom."""
    n_atoms = len(structure.atoms)
    component = [-1] * n_atoms
    molecules = []
    for start in range(n_atoms):
        if component[start] >= 0:
            continue
        component[start] = len(molecules)
        members, queue = [], deque([structure.atoms[start]])
        while queue:
            atom = queue.popleft()
            members.append(atom.idx)
            for partner in atom.bond_partners:
                if component[partner.idx] < 0:
                    component[partner.idx] = len(molecules)
                    queue.append(partner)
        molecules.append(sorted(members))

    split = []
    for members in molecules:
        by_chain = defaultdict(list)
        for index in members:
            residue = structure.atoms[index].residue
            by_chain[getattr(residue, "chain", "") or ""].append(index)
        split.extend(by_chain.values())
    return sorted(split, key=lambda members: members[0])



def _rigid_water(bonds: list) -> np.ndarray | None:
    """
    Returns the reference positions of a rigid three-site water in the
    atom order of its molecule, or None if the bonds are not an O-H-H
    triangle.
    """
    if len(bonds) != 3:
        return None
    atoms = sorted(
        {
            atom.idx: atom
            for bond in bonds
            for atom in (bond.atom1, bond.atom2)
        }.items()
    )
    if sorted(atom.atomic_number for _, atom in atoms) != [1, 1, 8]:
        return None

    local = {index: k for k, (index, _) in enumerate(atoms)}
    distance = {}
    for bond in bonds:
        if bond.type is None:
            raise Parm2PQError(
                f"The bond {bond.atom1.idx + 1}-{bond.atom2.idx + 1} has no "
                "parameters."
            )
        pair = _pair(local[bond.atom1.idx], local[bond.atom2.idx])
        distance[pair] = bond.type.req
    if len(distance) != 3:
        return None

    d_01, d_02, d_12 = distance[(0, 1)], distance[(0, 2)], distance[(1, 2)]
    x = (d_01**2 + d_02**2 - d_12**2) / (2 * d_01)
    y_squared = d_02**2 - x**2
    if y_squared <= 0.0:
        raise Parm2PQError(
            f"The water bonds {d_01}, {d_02} and {d_12} Angstrom of atoms "
            f"{atoms[0][0] + 1}-{atoms[-1][0] + 1} do not form a triangle."
        )
    return np.array(
        [[0.0, 0.0, 0.0], [d_01, 0.0, 0.0], [x, math.sqrt(y_squared), 0.0]]
    )



def _element(atom: Any) -> str:
    from parmed.periodic_table import Element  # pylint: disable=import-outside-toplevel
    return Element[atom.atomic_number]



class _Converter:  # pylint: disable=too-many-instance-attributes

    """Builds a :py:class:`PQForceFieldSetup` from a ParmEd structure."""

    def __init__(self, structure: Any, shake_h: bool):
        self.structure = structure
        self.shake_h = shake_h
        self.lj = _LennardJones(structure)
        if structure.adjusts:
            self.scale_14 = _scale_from_adjusts(structure, self.lj.rule)
        else:
            scale_coulomb, scale_vdw = _scale_from_dihedrals(
                structure, self.lj
            )
            self.scale_14 = (
                scale_coulomb, scale_vdw * lj_14_factor(structure)
            )
        self.excluded, self.pairs_14 = _nonbonded_pairs(structure)

        self.molecules = _molecules(structure)
        self.order = np.array([i for mol in self.molecules for i in mol])
        self.new_index = np.empty(len(self.order), dtype=int)
        self.new_index[self.order] = np.arange(len(self.order))
        self.molecule_of = np.empty(len(self.order), dtype=int)
        for index, members in enumerate(self.molecules):
            self.molecule_of[members] = index

        self.rigid = self._rigid_waters() if shake_h else {}
        self.mshake = {}

        self.intra = defaultdict(list)
        for pairs, value in ((self.excluded, 0), (self.pairs_14, -1)):
            for index_1, index_2 in pairs:
                molecule = self.molecule_of[index_1]
                if molecule == self.molecule_of[index_2]:
                    self.intra[molecule].append((index_1, index_2, value))

    def _rigid_waters(self) -> dict[int, np.ndarray]:
        """Returns the M-SHAKE reference of each rigid water molecule."""
        bonds_of = defaultdict(list)
        for bond in self.structure.bonds:
            bonds_of[self.molecule_of[bond.atom1.idx]].append(bond)
        rigid = {}
        for molecule, members in enumerate(self.molecules):
            if len(members) == 3:
                reference = _rigid_water(bonds_of[molecule])
                if reference is not None:
                    rigid[molecule] = reference
        return rigid

    def _intra_matrix(self, molecule: int) -> np.ndarray:
        members = self.molecules[molecule]
        n_atoms = len(members)
        local = {index: k for k, index in enumerate(members)}
        matrix = np.ones((n_atoms, n_atoms), dtype=np.int8)
        for index_1, index_2, value in self.intra[molecule]:
            k_1, k_2 = local[index_1], local[index_2]
            matrix[k_1, k_2] = matrix[k_2, k_1] = value
        return matrix

    def moltypes(self) -> tuple[list[PQMolType], np.ndarray]:
        """Returns the moltypes and the one-based moltype of each atom."""
        atoms = self.structure.atoms
        moltypes, index_of = [], {}
        atom_moltypes = np.empty(len(atoms), dtype=int)
        for molecule, members in enumerate(self.molecules):
            matrix = self._intra_matrix(molecule)
            per_atom = tuple(
                (
                    atoms[i].name,
                    str(atoms[i].type),
                    round(atoms[i].charge, 8),
                    self.lj.atom_vdw_types[i],
                    atoms[i].atomic_number,
                ) for i in members
            )
            key = (per_atom, matrix.tobytes())
            if key not in index_of:
                index_of[key] = len(moltypes) + 1
                moltypes.append(self._moltype(members, matrix))
            atom_moltypes[self.new_index[members]] = index_of[key]
            self._add_mshake(molecule, index_of[key])

        names = defaultdict(int)
        for moltype in moltypes:
            names[moltype.name] += 1
        counter = defaultdict(int)
        for moltype in moltypes:
            if names[moltype.name] > 1:
                counter[moltype.name] += 1
                moltype.name = f"{moltype.name}_{counter[moltype.name]}"
        return moltypes, atom_moltypes

    def _add_mshake(self, molecule: int, moltype: int) -> None:
        reference = self.rigid.get(molecule)
        if moltype not in self.mshake:
            self.mshake[moltype] = reference
            return
        known = self.mshake[moltype]
        if (reference is None) != (known is None) or (
            known is not None and not np.allclose(reference, known)
        ):
            raise Parm2PQError(
                f"The waters of moltype {moltype} are not all rigid with the "
                "same bond lengths; M-SHAKE needs one reference geometry per "
                "moltype."
            )

    def _moltype(self, members: list[int], matrix: np.ndarray) -> PQMolType:
        atoms = [self.structure.atoms[i] for i in members]
        residues = {atom.residue.idx for atom in atoms}
        name = atoms[0].residue.name if len(residues) == 1 else "MOL"
        name = "".join(name.split()) or "MOL"
        if name.lstrip("-").isdigit():
            name = f"MOL{name}"

        type_index = {}
        atom_types = [
            type_index.setdefault(str(atom.type), len(type_index))
            for atom in atoms
        ]
        intra_pairs = []
        for k in range(len(members)):
            row = matrix[k, k + 1:]
            partners = np.nonzero(row)[0]
            intra_pairs.append(
                [int(row[p]) * int(k + p + 2) for p in partners]
            )
        return PQMolType(
            name=name,
            elements=[_element(atom) for atom in atoms],
            atom_types=atom_types,
            charges=[float(atom.charge) for atom in atoms],
            vdw_types=[self.lj.atom_vdw_types[i] for i in members],
            intra_pairs=intra_pairs,
        )

    def bonded(self, setup: PQForceFieldSetup) -> None:
        """Fills the bonded types and terms of the setup."""
        builder = _BondedBuilder(self, setup)
        builder.bonds(self.structure.bonds, self.shake_h)
        builder.angles(self.structure.angles)
        builder.dihedrals(self.structure.dihedrals)
        builder.check_linkers()



class _BondedBuilder:

    """Bonded types and terms with the PQ linker marks."""

    def __init__(self, converter: _Converter, setup: PQForceFieldSetup):
        self.setup = setup
        self.new_index = converter.new_index
        molecule_of = converter.molecule_of
        self.rigid_atoms = {
            index
            for molecule in converter.rigid
            for index in converter.molecules[molecule]
        }
        self.cross_full = {
            pair
            for pair in converter.excluded
            if molecule_of[pair[0]] != molecule_of[pair[1]]
        }
        self.cross_14 = {
            pair
            for pair in converter.pairs_14
            if molecule_of[pair[0]] != molecule_of[pair[1]]
        }
        self.types = {}

    def _new(self, indices: list[int]) -> list[int]:
        return [int(self.new_index[i]) for i in indices]

    @staticmethod
    def _linker(indices: list[int], pairs: set) -> bool:
        """Marks the first term that connects a pair between molecules."""
        pair = _pair(indices[0], indices[-1])
        if pair in pairs:
            pairs.discard(pair)
            return True
        return False

    def _type_index(self, table: list, key: tuple) -> int:
        key = tuple(round(value, 10) for value in key)
        store = self.types.setdefault(id(table), {})
        if key not in store:
            store[key] = len(table)
            table.append(key)
        return store[key]

    def _is_rigid(self, indices: list[int]) -> bool:
        return all(index in self.rigid_atoms for index in indices)

    def bonds(self, bonds: Any, shake_h: bool) -> None:
        """
        Adds the bonds, or SHAKE constraints for bonds to hydrogen; M-SHAKE
        replaces the bonds of rigid water.
        """
        for bond in bonds:
            indices = [bond.atom1.idx, bond.atom2.idx]
            if self._is_rigid(indices):
                continue
            if bond.type is None:
                raise Parm2PQError(
                    f"The bond {indices[0] + 1}-{indices[1] + 1} has no "
                    "parameters."
                )
            is_h = 1 in (bond.atom1.atomic_number, bond.atom2.atomic_number)
            if shake_h and is_h:
                self.setup.shake.append((*self._new(indices), bond.type.req))
                continue
            index = self._type_index(
                self.setup.bond_types, (bond.type.req, 2 * bond.type.k)
            )
            linker = self._linker(indices, self.cross_full)
            self.setup.bonds.append((*self._new(indices), index, linker))

    def angles(self, angles: Any) -> None:
        """Adds the angles; M-SHAKE replaces the angle of rigid water."""
        for angle in angles:
            indices = [atom.idx for atom in _atoms_of(angle)]
            if self._is_rigid(indices):
                continue
            if angle.type is None:
                raise Parm2PQError(
                    f"The angle {'-'.join(str(i + 1) for i in indices)} has "
                    "no parameters."
                )
            index = self._type_index(
                self.setup.angle_types, (angle.type.theteq, 2 * angle.type.k)
            )
            linker = self._linker(indices, self.cross_full)
            self.setup.angles.append((*self._new(indices), index, linker))

    def dihedrals(self, dihedrals: Any) -> None:
        """Adds one proper or improper term per Fourier term."""
        for dihedral in dihedrals:
            indices = [atom.idx for atom in _atoms_of(dihedral)]
            new = self._new(indices)
            if dihedral.improper:
                table = self.setup.improper_types
            else:
                table = self.setup.dihedral_types
            for term in _dihedral_terms(dihedral):
                phase = (-term.phase) % 360.0
                index = self._type_index(table, (term.phi_k, term.per, phase))
                if dihedral.improper:
                    self.setup.impropers.append((*new, index))
                else:
                    linker = self._linker(indices, self.cross_14)
                    self.setup.dihedrals.append((*new, index, linker))

    def check_linkers(self) -> None:
        """Raises if a pair between molecules has no carrier term."""
        missing = sorted(self.cross_full | self.cross_14)
        if missing:
            index_1, index_2 = missing[0]
            raise Parm2PQError(
                f"{len(missing)} excluded or 1-4 atom pairs between molecules "
                "have no bond, angle or proper dihedral to carry the PQ "
                f"linker correction (first: atoms {index_1 + 1} and "
                f"{index_2 + 1})."
            )



def _box(structure: Any, positions: np.ndarray, rcoulomb: Real) -> np.ndarray:
    if structure.box is not None:
        return np.asarray(structure.box, dtype=float)
    extent = positions.max(axis=0) - positions.min(axis=0)
    edge = float(extent.max()) + 2.0 * rcoulomb
    return np.array([edge, edge, edge, 90.0, 90.0, 90.0])



def convert_structure(
    structure: Any,
    shake_h: bool = False,
    rcoulomb: Real = 12.0,
    source: str = "",
) -> PQForceFieldSetup:
    """
    Converts a parametrized ParmEd structure into a PQ setup.

    Parameters
    ----------
    structure : parmed.Structure
        The parametrized structure with coordinates.
    shake_h : bool, optional
        Constrain bonds to hydrogen with SHAKE instead of a bond term, by
        default False. Rigid three-site water (with an H-H bond) is kept
        rigid with M-SHAKE instead.
    rcoulomb : float, optional
        The Coulomb cutoff in Angstrom, used for the vacuum box of a
        structure without box, by default 12.0.
    source : str, optional
        A description of the source for the file headers.

    Returns
    -------
    PQForceFieldSetup
        The PQ setup.

    Raises
    ------
    UnsupportedForceFieldTermError
        If the force field contains a term PQ cannot evaluate.
    Parm2PQError
        If the structure has no coordinates or misses parameters.
    """
    check_supported(structure)
    if structure.coordinates is None:
        raise Parm2PQError(
            "The structure has no coordinates; give a coordinate file."
        )

    converter = _Converter(structure, shake_h)
    moltypes, atom_moltypes = converter.moltypes()
    order = converter.order
    positions = np.asarray(structure.coordinates, dtype=float)[order]

    velocities = getattr(structure, "velocities", None)
    if velocities is not None:
        velocities = np.asarray(velocities, dtype=float)[order] * 1e12

    setup = PQForceFieldSetup(
        moltypes=moltypes,
        atom_moltypes=atom_moltypes,
        elements=[_element(structure.atoms[i]) for i in order],
        positions=positions,
        velocities=velocities,
        box=_box(structure, positions, rcoulomb),
        atom_order=order,
        scale_14_coulomb=converter.scale_14[0],
        scale_14_vdw=converter.scale_14[1],
        lj_pairs=converter.lj.all_pairs(),
        mshake={
            moltype: reference
            for moltype, reference in converter.mshake.items()
            if reference is not None
        },
        source=source,
        periodic=structure.box is not None,
    )
    converter.bonded(setup)
    return setup
