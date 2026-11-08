"""
A module containing the centers of mass of molecules for the RDF analysis.

The centers are appended to the positions of every frame as additional rows,
so that the distance kernels of the RDF analysis (compiled and fallback)
work on them like on atoms.
"""

import numpy as np

from PQAnalysis.atomic_system._center_of_mass import center_of_mass_groups
from PQAnalysis.core import Cell



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
