"""
A module containing the centers of mass of molecules for the analyses.

The centers are appended to the positions of every frame as additional rows,
so that the kernels of the RDF and MSD analyses (compiled and fallback) work
on them like on atoms.
"""

import numpy as np

from PQAnalysis.atomic_system._center_of_mass import (
    center_of_mass_groups,
    placeholder_atom_names,
)
from PQAnalysis.core import Cell
from PQAnalysis.topology import Topology



class MoleculeCenters:

    """
    The centers of mass of molecules that are appended to the atoms.

    The positions of a frame with ``n_atoms`` atoms are extended by one row
    per molecule: the center of mass of the atoms ``atom_indices`` of this
    molecule. The rows are in ascending order of the molecule ids, row
    ``n_atoms + k`` belongs to the k-th smallest molecule id.
    """

    def __init__(
        self,
        n_atoms: int,
        atom_indices: np.ndarray,
        molecule_ids: np.ndarray,
        masses: np.ndarray,
    ) -> None:
        """
        Parameters
        ----------
        n_atoms : int
            The number of atoms of every frame.
        atom_indices : np.ndarray
            The indices of the atoms that define the centers.
        molecule_ids : np.ndarray
            The molecule of each of these atoms.
        masses : np.ndarray
            The mass of each of these atoms.
        """

        self.n_atoms = n_atoms
        self.atom_indices = np.asarray(atom_indices, dtype=np.intp)
        self.molecule_ids = np.asarray(molecule_ids)
        self.masses = np.asarray(masses, dtype=float)
        self.molecules = np.unique(self.molecule_ids)

    @classmethod
    def from_topology(
        cls,
        topology: Topology,
        atoms: np.ndarray,
        logger,
        exception,
    ) -> "MoleculeCenters":
        """
        Creates the centers of the molecules of the given atoms.

        The molecules are the residues of the topology and the atoms define
        the center of their molecule.

        Parameters
        ----------
        topology : Topology
            The topology with the residues.
        atoms : np.ndarray
            The indices of the atoms that define the centers.
        logger : logging.Logger
            The logger of the analysis, used for the errors and the warning.
        exception : type
            The exception class of the analysis that is raised for errors.

        Returns
        -------
        MoleculeCenters
            The centers, one per molecule that has a given atom.

        Raises
        ------
        exception
            If the topology has no residues or if an atom has no mass.
        """

        if not topology.residues:
            logger.error(
                (
                    "The centers of mass need the molecules of the system, "
                    "but the topology has no residues. Provide a restart "
                    "file together with a moldescriptor file."
                ),
                exception=exception
            )

        topology_atoms = topology.atoms

        if any(topology_atoms[index].mass is None for index in atoms):
            logger.error(
                (
                    "The center of mass needs the mass of every selected "
                    "atom, but some selected atoms have no mass."
                ),
                exception=exception
            )

        placeholders = placeholder_atom_names(topology_atoms, atoms)
        if placeholders:
            logger.warning(
                (
                    "The selection contains atoms with a placeholder mass "
                    f"({', '.join(placeholders)}). They dominate the center "
                    "of mass of their molecule; exclude them with the "
                    "selection."
                )
            )

        return cls(
            n_atoms=topology.n_atoms,
            atom_indices=atoms,
            molecule_ids=topology.residue_numbers[atoms],
            masses=np.array([topology_atoms[index].mass for index in atoms]),
        )

    @property
    def n_centers(self) -> int:
        """int: The number of centers of mass."""
        return len(self.molecules)

    def rows(self, molecule_ids: np.ndarray) -> np.ndarray:
        """
        Returns the rows of the extended positions that hold the centers
        of the given molecules.

        Parameters
        ----------
        molecule_ids : np.ndarray
            The molecules, which all have to be molecules of the centers.

        Returns
        -------
        np.ndarray
            The row indices into the extended positions.
        """

        return self.n_atoms + np.searchsorted(self.molecules, molecule_ids)

    def extend(self, positions: np.ndarray, cell: Cell) -> np.ndarray:
        """
        Appends the centers of mass to the positions of a frame.

        Parameters
        ----------
        positions : np.ndarray
            The positions of the atoms with shape (n_atoms, 3).
        cell : Cell
            The cell of the frame.

        Returns
        -------
        np.ndarray
            A new C-contiguous float64 array with shape
            (n_atoms + n_centers, 3).

        Raises
        ------
        AtomicSystemError
            If a molecule spans half of the box or more.
        """

        extended = np.empty((self.n_atoms + self.n_centers, 3))
        extended[:self.n_atoms] = positions
        extended[self.n_atoms:] = center_of_mass_groups(
            extended[self.atom_indices],
            self.masses,
            self.molecule_ids,
            cell,
            check_extent=True,
        )

        return extended
