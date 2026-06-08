"""Pure Python fallback for trajectory line parsing."""

import numpy as np


def _split_xyz_line(line: str) -> tuple[str, float, float, float]:
    fields = line.split()
    if len(fields) < 4:
        raise ValueError("Could not parse line")

    try:
        return fields[0], float(fields[1]), float(fields[2]), float(fields[3])
    except ValueError as exc:
        raise ValueError("Could not parse line") from exc


def process_lines_with_atoms(input, n_atoms: int):
    """
    Parse atom labels and xyz coordinates from trajectory lines.
    """

    xyz = np.zeros((n_atoms, 3), dtype=np.float32)
    atoms = [None] * n_atoms

    for i in range(n_atoms):
        atom, x, y, z = _split_xyz_line(input[i])
        xyz[i, 0] = x
        xyz[i, 1] = y
        xyz[i, 2] = z
        atoms[i] = atom

    return atoms, xyz


def process_lines(input, n_atoms: int):
    """
    Parse xyz coordinates from trajectory lines.
    """

    xyz = np.zeros((n_atoms, 3), dtype=np.float32)

    for i in range(n_atoms):
        _, x, y, z = _split_xyz_line(input[i])
        xyz[i, 0] = x
        xyz[i, 1] = y
        xyz[i, 2] = z

    return xyz
