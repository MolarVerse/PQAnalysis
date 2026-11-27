"""
A module containing a periodic cell list for fast overlap checks.
"""

import numpy as np

from PQAnalysis.core import Cell
from PQAnalysis.types import Np2DNumberArray, PositiveReal

# upper bound for the number of grid cells to keep memory bounded
_MAX_GRID_CELLS = 2**21

# keeps FP noise in the binning from splitting a cutoff sphere over two cells
_CELL_SIZE_MARGIN = 1.0 + 1e-6

_NEIGHBOUR_OFFSETS = np.array(
    np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij")
).reshape(3, -1).T



class OverlapGrid:  # pylint: disable=too-many-instance-attributes

    """
    A periodic cell list that answers whether new positions overlap
    with the stored positions.

    The decision is the same as checking all pairs with
    :py:func:`~PQAnalysis.core.api.distance`: positions overlap if any
    minimum image distance is smaller than or equal to the cutoff. Each
    grid cell is at least one cutoff wide in every direction, so only
    the 27 neighbouring cells have to be searched. Positions can be
    added one molecule at a time.
    """

    def __init__(
        self,
        cell: Cell,
        cutoff: PositiveReal,
        positions: Np2DNumberArray | None = None,
    ) -> None:
        """
        Parameters
        ----------
        cell : Cell
            The periodic cell. Must not be a vacuum cell.
        cutoff : PositiveReal
            The distance cutoff.
        positions : Np2DNumberArray | None, optional
            Initial positions, by default None
        """
        self._cutoff = cutoff
        self._box_lengths = np.asarray(cell.box_lengths, dtype=float)
        self._box_matrix = cell.box_matrix
        self._inverse_box_matrix = cell.inverse_box_matrix
        self._orthorhombic = bool(
            cell.alpha == 90 and cell.beta == 90 and cell.gamma == 90
        )

        self._n_cells = self._grid_shape(cell, cutoff)
        n_total = int(np.prod(self._n_cells))
        self._strides = np.array(
            [self._n_cells[1] * self._n_cells[2], self._n_cells[2], 1]
        )
        self._to_grid = self._inverse_box_matrix.T * self._n_cells

        self._table = np.full((n_total, 4), -1, dtype=np.int64)
        self._counts = np.zeros(n_total, dtype=np.int64)
        # the last row stays NaN so that empty slots (-1) never overlap
        self._pos = np.full((1, 3), np.nan)
        self._n_pos = 0

        if positions is not None:
            self.add(positions)

    @staticmethod
    def _grid_shape(cell: Cell, cutoff: PositiveReal) -> np.ndarray:
        box = cell.box_matrix
        a, b, c = box[:, 0], box[:, 1], box[:, 2]
        volume = abs(np.dot(a, np.cross(b, c)))
        widths = volume / np.array(
            [
                np.linalg.norm(np.cross(b, c)),
                np.linalg.norm(np.cross(c, a)),
                np.linalg.norm(np.cross(a, b)),
            ]
        )

        n_cells = np.floor(widths / (cutoff * _CELL_SIZE_MARGIN))
        n_cells = np.clip(n_cells, 1, None)

        total = np.prod(n_cells)
        if total > _MAX_GRID_CELLS:
            n_cells = np.floor(n_cells / np.cbrt(total / _MAX_GRID_CELLS))
            n_cells = np.clip(n_cells, 1, None)

        return n_cells.astype(np.int64)

    @property
    def n_positions(self) -> int:
        """int: The number of stored positions."""
        return self._n_pos

    @property
    def positions(self) -> Np2DNumberArray:
        """Np2DNumberArray: The stored positions."""
        return self._pos[:self._n_pos]

    def _cell_indices(self, positions: Np2DNumberArray) -> np.ndarray:
        return np.floor(positions @ self._to_grid).astype(np.int64)

    def _flat(self, indices: np.ndarray) -> np.ndarray:
        return (indices % self._n_cells) @ self._strides

    def add(self, positions: Np2DNumberArray) -> None:
        """
        Store new positions.

        Parameters
        ----------
        positions : Np2DNumberArray
            The positions to add.
        """
        positions = np.reshape(np.asarray(positions, dtype=float), (-1, 3))
        n_new = len(positions)
        if n_new == 0:
            return

        if self._n_pos + n_new + 1 > len(self._pos):
            size = max(2 * len(self._pos), self._n_pos + n_new + 1)
            grown = np.full((size, 3), np.nan)
            grown[:self._n_pos] = self._pos[:self._n_pos]
            self._pos = grown

        self._pos[self._n_pos:self._n_pos + n_new] = positions

        flat = self._flat(self._cell_indices(positions))
        for offset, cell_id in enumerate(flat):
            count = self._counts[cell_id]
            if count == self._table.shape[1]:
                self._grow_table()
            self._table[cell_id, count] = self._n_pos + offset
            self._counts[cell_id] = count + 1

        self._n_pos += n_new

    def _grow_table(self) -> None:
        n_cells, capacity = self._table.shape
        grown = np.full((n_cells, 2 * capacity), -1, dtype=np.int64)
        grown[:, :capacity] = self._table
        self._table = grown

    def _image(self, delta: np.ndarray) -> np.ndarray:
        # same arithmetic as Cell.image so decisions match distance()
        if self._orthorhombic:
            return delta - self._box_lengths * np.round(
                delta / self._box_lengths
            )

        fractional = delta @ self._inverse_box_matrix.T
        fractional -= np.round(fractional)
        return fractional @ self._box_matrix.T

    def overlaps(self, positions: Np2DNumberArray) -> bool:
        """
        Check if any of the given positions is within the cutoff of a
        stored position.

        Parameters
        ----------
        positions : Np2DNumberArray
            The positions to check.

        Returns
        -------
        bool
            True if any minimum image distance is smaller than or
            equal to the cutoff.
        """
        if self._n_pos == 0:
            return False

        positions = np.reshape(positions, (-1, 3))
        neighbours = self._flat(
            self._cell_indices(positions)[:, None, :] + _NEIGHBOUR_OFFSETS
        )
        candidates = self._table[neighbours].reshape(len(positions), -1)

        delta = positions[:, None, :] - self._pos[candidates]
        distances = np.linalg.norm(self._image(delta), axis=-1)

        return bool(np.any(distances <= self._cutoff))
