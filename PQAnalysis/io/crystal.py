"""Prepare ordered periodic crystal geometry from CIF."""

from math import prod
import warnings

import numpy as np

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell



def prepare_crystal(
    filename: str,
    repeats: tuple[int, int, int] = (1, 1, 1),
    *,
    block: str | None = None,
    max_atoms: int = 100_000,
) -> AtomicSystem:
    """Read an ordered CIF and build a centered periodic supercell.

    Parameters
    ----------
    filename : str
        Input CIF path.
    repeats : tuple of int
        Positive replication counts along the three lattice vectors.
    block : str or None
        CIF data block name without the ``data_`` prefix. Required when there
        is more than one structure block.
    max_atoms : int
        Maximum atom count after symmetry expansion and replication.

    Returns
    -------
    AtomicSystem
        Geometry only, in angstroms and in the PQ cell convention.

    Raises
    ------
    ValueError
        If block selection is ambiguous, the geometry is invalid, sites
        are disordered or coincident, or the atom budget is exceeded.
    ImportError
        If ASE is not installed.
    """
    from ase.geometry import cellpar_to_cell  # pylint: disable=import-outside-toplevel

    if len(repeats) != 3 or any(
        isinstance(n, bool) or not isinstance(n, int) or n < 1 for n in repeats
    ):
        raise ValueError("Specify three positive integer repeats.")
    if max_atoms < 1:
        raise ValueError("The atom limit must be positive.")

    crystal = _read_ordered_cif(filename, block)
    if len(crystal) * prod(repeats) > max_atoms:
        raise ValueError("The requested supercell exceeds the atom limit.")

    supercell = crystal.repeat(repeats)
    parameters = supercell.cell.cellpar()
    fractional = supercell.get_scaled_positions(wrap=False)
    fractional -= np.floor(fractional + 0.5)
    return AtomicSystem(
        atoms=[Atom(symbol) for symbol in supercell.get_chemical_symbols()],
        pos=fractional @ cellpar_to_cell(parameters),
        cell=Cell(*(float(value) for value in parameters)),
    )



def _read_ordered_cif(filename, block):
    """Read one ordered, fully periodic CIF structure."""
    # ASE is optional for unrelated workflows.
    # pylint: disable=import-outside-toplevel
    from ase.io.cif import parse_cif
    from ase.neighborlist import neighbor_list

    candidates = [
        entry for entry in parse_cif(filename)
        if entry.has_structure() and (block is None or entry.name == block)
    ]
    if len(candidates) != 1:
        raise ValueError("Select exactly one CIF structure block.")
    entry = candidates[0]
    if "_atom_site_occupancy" in entry and not np.allclose(
        np.asarray(entry["_atom_site_occupancy"], dtype=float),
        1.0,
        rtol=0.0,
        atol=1e-8,
    ):
        raise ValueError("Resolve partial occupancies before building a cell.")

    cell = entry.get_cell()
    if not np.isfinite(cell.array).all() or cell.volume <= 0:
        raise ValueError("The CIF must describe a finite, nonsingular cell.")

    # ASE otherwise warns and can collapse different source sites.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "error", message="scaled_positions.*equivalent"
        )
        try:
            crystal = entry.get_atoms(fractional_occupancies=False)
        except UserWarning as error:
            raise ValueError(str(error)) from error

    if len(crystal) == 0 or not np.isfinite(crystal.positions).all():
        raise ValueError("The CIF must contain finite atomic positions.")
    if neighbor_list("i", crystal, 0.01).size:
        raise ValueError("The CIF contains coincident periodic sites.")
    return crystal
