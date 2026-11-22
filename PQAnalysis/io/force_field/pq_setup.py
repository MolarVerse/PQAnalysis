"""
The PQ force-field setup produced by :py:func:`parm2pq` and its writers.

All atom indices stored here are zero-based positions in the PQ atom
order; the writers convert them to the one-based indices of the PQ
setup files.
"""

from dataclasses import dataclass, field
from numbers import Real
from pathlib import Path

import numpy as np

#: default file names of the written PQ setup files
DEFAULT_FILE_NAMES = {
    "moldescriptor": "moldescriptor.dat",
    "topology": "topology.top",
    "parameter": "parameter.prm",
    "intra_nonbonded": "intra_nonbonded.dat",
    "start": "start.rst",
    "input": "pq_input.in",
    "mshake": "mshake.dat",
}



def _fmt(value: Real) -> str:
    """Formats a float with 12 significant digits and without ``-0``."""
    text = f"{value:.12g}"
    return "0" if text == "-0" else text



@dataclass
class PQMolType:

    """
    One PQ moltype: the per-atom data of the moldescriptor file and the
    intramolecular non-bonded pairs.
    """

    #: the moltype name
    name: str
    #: the element symbol of each atom; PQ derives the masses from it
    elements: list[str]
    #: the atom type index of each atom within the moltype
    atom_types: list[int]
    #: the partial charge of each atom in e
    charges: list[float]
    #: the global van der Waals index of each atom
    vdw_types: list[int]
    #: for each atom (zero-based) the interacting partners (one-based) with
    #: a larger index; negative partners are scaled 1-4 pairs
    intra_pairs: list[list[int]]

    @property
    def n_atoms(self) -> int:
        """int: The number of atoms of the moltype."""
        return len(self.elements)



@dataclass
class PQForceFieldSetup:  # pylint: disable=too-many-instance-attributes

    """
    A complete PQ force-field setup. Atom indices are zero-based positions
    in the PQ atom order; the writers add one.
    """

    #: the moltypes; moltype ``k`` in PQ is ``moltypes[k - 1]``
    moltypes: list[PQMolType]
    #: the one-based moltype of each atom
    atom_moltypes: np.ndarray
    #: the element symbol of each atom
    elements: list[str]
    #: the positions in Angstrom, shape (n_atoms, 3)
    positions: np.ndarray
    #: the velocities in Angstrom/s, shape (n_atoms, 3), or None
    velocities: np.ndarray | None
    #: the box lengths in Angstrom and angles in degrees
    box: np.ndarray
    #: the index in the source structure of each PQ atom
    atom_order: np.ndarray
    #: the 1-4 Coulomb scaling factor
    scale_14_coulomb: Real
    #: the 1-4 van der Waals scaling factor
    scale_14_vdw: Real
    #: equilibrium distance in Angstrom and PQ force constant (1/2 k)
    bond_types: list = field(default_factory=list)
    #: equilibrium angle in degrees and PQ force constant (1/2 k)
    angle_types: list = field(default_factory=list)
    #: force constant, periodicity and PQ phase in degrees
    dihedral_types: list = field(default_factory=list)
    #: force constant, periodicity and PQ phase in degrees
    improper_types: list = field(default_factory=list)
    #: global van der Waals indices, c6 and c12
    lj_pairs: list = field(default_factory=list)
    #: atom indices, type index and linker flag
    bonds: list = field(default_factory=list)
    #: atom indices, type index and linker flag
    angles: list = field(default_factory=list)
    #: atom indices, type index and linker flag
    dihedrals: list = field(default_factory=list)
    #: atom indices and type index
    impropers: list = field(default_factory=list)
    #: constrained atom pair and distance in Angstrom
    shake: list = field(default_factory=list)
    #: M-SHAKE reference positions in Angstrom by one-based moltype
    mshake: dict = field(default_factory=dict)
    #: a description of the source force field for the file headers
    source: str = ""

    @property
    def n_atoms(self) -> int:
        """int: The number of atoms."""
        return len(self.elements)

    def _header(self) -> list[str]:
        lines = ["# written by PQAnalysis parm2pq"]
        if self.source:
            lines.append(f"# source: {self.source}")
        return lines

    def moldescriptor_lines(self) -> list[str]:
        """
        Returns the lines of the moldescriptor file.

        Returns
        -------
        list[str]
            The lines without line breaks.
        """
        lines = self._header()
        for index, moltype in enumerate(self.moltypes, start=1):
            charge = sum(moltype.charges)
            lines.append("")
            lines.append(f"# MOLTYPE {index}")
            lines.append(
                f"{moltype.name} {moltype.n_atoms} {_fmt(round(charge, 8))}"
            )
            for element, atom_type, charge, vdw_type in zip(
                moltype.elements,
                moltype.atom_types,
                moltype.charges,
                moltype.vdw_types,
            ):
                lines.append(
                    f"  {element:<3} {atom_type:>4} {_fmt(charge):>16} "
                    f"{vdw_type:>4}"
                )
        return lines

    def intra_nonbonded_lines(self) -> list[str]:
        """
        Returns the lines of the intra non-bonded file.

        Returns
        -------
        list[str]
            The lines without line breaks.
        """
        lines = self._header()
        for index, moltype in enumerate(self.moltypes, start=1):
            lines.append(f"# {moltype.name}")
            lines.append(f"{index}")
            for atom, partners in enumerate(moltype.intra_pairs, start=1):
                if partners:
                    lines.append(
                        f"{atom} " + " ".join(str(p) for p in partners)
                    )
            lines.append("END")
        return lines

    def topology_lines(self) -> list[str]:
        """
        Returns the lines of the topology file.

        Returns
        -------
        list[str]
            The lines without line breaks.
        """
        lines = self._header()

        def section(name, entries, formatter):
            if not entries:
                return
            lines.append("")
            lines.append(name)
            lines.extend(formatter(entry) for entry in entries)
            lines.append("END")

        def linked(atoms, type_index, linker):
            text = " ".join(str(atom + 1) for atom in atoms)
            text += f" {type_index}"
            return text + " *" if linker else text

        section(
            "SHAKE",
            self.shake, lambda e: f"{e[0] + 1} {e[1] + 1} "
            f"{_fmt(e[2])}"
        )
        section("BONDS", self.bonds, lambda e: linked(e[:2], e[2], e[3]))
        section("ANGLES", self.angles, lambda e: linked(e[:3], e[3], e[4]))
        section(
            "DIHEDRALS",
            self.dihedrals,
            lambda e: linked(e[:4], e[4], e[5]),
        )
        section(
            "IMPROPERS",
            self.impropers,
            lambda e: linked(e[:4], e[4], False),
        )
        return lines

    def mshake_lines(self) -> list[str]:
        """
        Returns the lines of the M-SHAKE file: one reference geometry per
        rigid moltype. PQ reads no comment or empty lines in this file.

        Returns
        -------
        list[str]
            The lines without line breaks.
        """
        lines = []
        for moltype, positions in sorted(self.mshake.items()):
            elements = self.moltypes[moltype - 1].elements
            lines.append(f"{len(elements)}")
            lines.append(f"moltype = {moltype};")
            for element, position in zip(elements, positions):
                lines.append(
                    f"{element:<3} " + " ".join(f"{x:.10f}" for x in position)
                )
        return lines

    def parameter_lines(self) -> list[str]:
        """
        Returns the lines of the parameter file.

        Returns
        -------
        list[str]
            The lines without line breaks.
        """
        lines = self._header()
        lines.append("#     dummy x5  scale_14_coulomb  scale_14_vdw")
        lines.append(
            f"TYPES 0 0 0 0 0 {_fmt(self.scale_14_coulomb)} "
            f"{_fmt(self.scale_14_vdw)}"
        )

        def section(name, entries):
            if not entries:
                return
            lines.append("")
            lines.append(name)
            for index, entry in enumerate(entries):
                values = " ".join(_fmt(value) for value in entry)
                lines.append(f"{index} {values}")
            lines.append("END")

        section("BONDS", self.bond_types)
        section("ANGLES", self.angle_types)
        section("DIHEDRALS", self.dihedral_types)
        section("IMPROPERS", self.improper_types)

        lines.append("")
        lines.append("NONCOULOMBICS LJ")
        for type_1, type_2, c6, c12 in self.lj_pairs:
            lines.append(f"{type_1} {type_2} {_fmt(c6)} {_fmt(c12)}")
        lines.append("END")
        return lines

    def restart_lines(self) -> list[str]:
        """
        Returns the lines of the PQ restart (start) file.

        Returns
        -------
        list[str]
            The lines without line breaks.
        """
        lines = ["Box " + " ".join(_fmt(value) for value in self.box)]
        velocities = self.velocities
        if velocities is None:
            velocities = np.zeros_like(self.positions)

        counter = 0
        previous = None
        for i in range(self.n_atoms):
            moltype = int(self.atom_moltypes[i])
            counter = counter + 1 if moltype == previous else 1
            previous = moltype
            if counter > self.moltypes[moltype - 1].n_atoms:
                counter = 1
            pos = " ".join(f"{x:.10f}" for x in self.positions[i])
            vel = " ".join(f"{v:.8e}" for v in velocities[i])
            lines.append(
                f"{self.elements[i]} {counter} {moltype} {pos} {vel} 0 0 0"
            )
        return lines

    def input_lines(
        self,
        file_names: dict[str, str],
        rcoulomb: Real,
    ) -> list[str]:
        """
        Returns the lines of the PQ input snippet.

        Parameters
        ----------
        file_names : dict[str, str]
            The names of the written setup files.
        rcoulomb : float
            The Coulomb cutoff in Angstrom.

        Returns
        -------
        list[str]
            The lines without line breaks.
        """
        lines = self._header()
        lines += [
            "# add the run keys (nstep, timestep, thermostat, ...)",
            "jobtype = mm-md;",
            "force-field = on;",
            f"rcoulomb = {_fmt(rcoulomb)};",
            "long_range = none;",
            f"start_file = {file_names['start']};",
            f"moldescriptor_file = {file_names['moldescriptor']};",
            f"topology_file = {file_names['topology']};",
            f"parameter_file = {file_names['parameter']};",
            f"intra-nonBonded_file = {file_names['intra_nonbonded']};",
        ]
        if self.mshake:
            lines.append("shake = mshake;")
            lines.append(f"mshake_file = {file_names['mshake']};")
        elif self.shake:
            lines.append("shake = on;")
        return lines

    def write(
        self,
        output_dir: str | Path,
        rcoulomb: Real,
        overwrite: bool = False,
    ) -> dict[str, Path]:
        """
        Writes all PQ setup files into a directory.

        Parameters
        ----------
        output_dir : str | Path
            The directory; it is created if needed.
        rcoulomb : float
            The Coulomb cutoff in Angstrom written to the input snippet.
        overwrite : bool, optional
            Whether existing files may be replaced, by default False.

        Returns
        -------
        dict[str, Path]
            The written files by kind.

        Raises
        ------
        FileExistsError
            If a file exists and overwrite is False.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = {
            kind: output_dir / name
            for kind, name in DEFAULT_FILE_NAMES.items()
            if kind != "mshake" or self.mshake
        }

        if not overwrite:
            existing = [str(path) for path in paths.values() if path.exists()]
            if existing:
                raise FileExistsError(
                    "parm2pq does not replace existing files without "
                    f"--overwrite: {', '.join(existing)}"
                )

        contents = {
            "moldescriptor": self.moldescriptor_lines(),
            "topology": self.topology_lines(),
            "parameter": self.parameter_lines(),
            "intra_nonbonded": self.intra_nonbonded_lines(),
            "start": self.restart_lines(),
            "input": self.input_lines(DEFAULT_FILE_NAMES, rcoulomb),
        }
        if self.mshake:
            contents["mshake"] = self.mshake_lines()
        for kind, lines in contents.items():
            paths[kind].write_text("\n".join(lines) + "\n", encoding="utf-8")
        return paths
