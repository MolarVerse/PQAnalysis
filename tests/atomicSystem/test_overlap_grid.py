"""
Tests for the OverlapGrid cell list against the all-pairs distance check.
"""

import numpy as np
import pytest

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.atomic_system._overlap_grid import OverlapGrid
from PQAnalysis.atomic_system import _overlap_grid
from PQAnalysis.core import Atom, Cell, distance

from . import pytestmark  # pylint: disable=unused-import

CELLS = [
    Cell(10.0, 10.0, 10.0),
    Cell(7.3, 11.1, 9.4),
    Cell(2.5, 2.5, 2.5),
    Cell(9.0, 8.0, 10.0, 80.0, 100.0, 70.0),
    Cell(6.0, 9.0, 7.0, 60.0, 75.0, 110.0),
]



def brute_force_overlaps(stored, query, cell, cutoff):
    if len(stored) == 0:
        return False

    return bool(np.any(distance(stored, query, cell) <= cutoff))



def random_positions(rng, cell, n):
    # spans two box lengths so positions outside the box are covered
    fractional = rng.uniform(-1.0, 1.0, size=(n, 3))
    return fractional @ cell.box_matrix.T



class TestOverlapGrid:

    @pytest.mark.parametrize("cell", CELLS)
    @pytest.mark.parametrize("cutoff", [0.7, 1.0, 1.9])
    def test_matches_brute_force(self, cell, cutoff):
        rng = np.random.default_rng(1)
        # density giving a mix of overlapping and free queries
        n_stored = max(2, round(0.06 * cell.volume / cutoff**3))
        stored = random_positions(rng, cell, n_stored)
        grid = OverlapGrid(cell, cutoff, stored)

        decisions = []
        for _ in range(300):
            query = random_positions(rng, cell, 3)
            expected = brute_force_overlaps(stored, query, cell, cutoff)
            assert grid.overlaps(query) == expected
            decisions.append(expected)

        assert any(decisions)
        if 2 * cutoff < min(cell.box_lengths):
            assert not all(decisions)

    @pytest.mark.parametrize("cell", CELLS)
    def test_incremental_add_matches_brute_force(self, cell):
        rng = np.random.default_rng(2)
        cutoff = 1.0
        grid = OverlapGrid(cell, cutoff)
        stored = np.zeros((0, 3))

        assert not grid.overlaps(random_positions(rng, cell, 3))

        for _ in range(40):
            new = random_positions(rng, cell, 3)
            grid.add(new)
            stored = np.concatenate([stored, new])

            query = random_positions(rng, cell, 3)
            assert grid.overlaps(query) == brute_force_overlaps(
                stored, query, cell, cutoff
            )

        assert grid.n_positions == len(stored)
        assert np.array_equal(grid.positions, stored)

    @pytest.mark.parametrize("cell", CELLS)
    def test_periodic_edges(self, cell):
        cutoff = 1.0
        corner = np.zeros(3) + 1e-9
        stored = np.array([corner @ cell.box_matrix.T])
        grid = OverlapGrid(cell, cutoff, stored)

        for axis in range(3):
            unit = cell.box_matrix[:, axis] / np.linalg.norm(
                cell.box_matrix[:, axis]
            )
            for shift in (0.5, 0.999, 1.001, 1.5):
                for sign in (-1, 1):
                    query = (stored[0] + sign * shift * unit)[None]
                    assert grid.overlaps(query) == brute_force_overlaps(
                        stored, query, cell, cutoff
                    )

                    shifted = query + cell.box_matrix[:, axis]
                    assert grid.overlaps(shifted) == brute_force_overlaps(
                        stored, shifted, cell, cutoff
                    )

    def test_exact_cutoff_counts_as_overlap(self):
        cell = Cell(10.0, 10.0, 10.0)
        grid = OverlapGrid(cell, 1.0, np.array([[4.5, 0.0, 0.0]]))

        assert grid.overlaps(np.array([[-4.5, 0.0, 0.0]]))
        assert not grid.overlaps(np.array([[-4.0, 0.0, 0.0]]))

    def test_dense_cells_grow_the_table(self):
        cell = Cell(10.0, 10.0, 10.0)
        stored = np.full((20, 3), 0.25) + np.arange(20)[:, None] * 1e-3
        grid = OverlapGrid(cell, 1.0, stored)

        assert grid.overlaps(np.array([[0.9, 0.25, 0.25]]))
        assert not grid.overlaps(np.array([[2.0, 0.25, 0.25]]))

    def test_grid_size_is_capped(self, monkeypatch):
        monkeypatch.setattr(_overlap_grid, "_MAX_GRID_CELLS", 64)
        cell = Cell(50.0, 50.0, 50.0)
        rng = np.random.default_rng(3)
        stored = random_positions(rng, cell, 200)
        grid = OverlapGrid(cell, 1.5, stored)

        assert np.prod(grid._n_cells) <= 64

        for _ in range(100):
            query = random_positions(rng, cell, 3)
            assert grid.overlaps(query) == brute_force_overlaps(
                stored, query, cell, 1.5
            )



class TestFitAtomicSystem:

    @staticmethod
    def water():
        return AtomicSystem(
            atoms=[Atom("O"), Atom("H"), Atom("H")],
            pos=np.array(
                [
                    [0.0, 0.0, 0.0], [0.7572, 0.5865, 0.0],
                    [-0.7572, 0.5865, 0.0]
                ]
            ),
        )

    @staticmethod
    def framework():
        rng = np.random.default_rng(4)
        cell = Cell(12.0, 12.0, 12.0)
        return AtomicSystem(
            atoms=[Atom("C")] * 30,
            pos=rng.uniform(-6.0, 6.0, size=(30, 3)),
            cell=cell,
        )

    def test_fitted_molecules_keep_the_cutoff(self, monkeypatch):
        monkeypatch.setenv("CONSTANT_SEED", "7")
        framework = self.framework()
        water = self.water()

        fitted = framework.fit_atomic_system(
            water, number_of_additions=20, distance_cutoff=1.2
        )

        assert len(fitted) == 20
        placed = [framework.pos]
        for system in fitted:
            assert np.allclose(
                distance(system.pos, system.pos, framework.cell),
                distance(water.pos, water.pos),
            )
            assert np.all(
                distance(np.concatenate(placed), system.pos, framework.cell) >
                1.2
            )
            placed.append(system.pos)

    def test_constant_seed_is_reproducible(self, monkeypatch):
        monkeypatch.setenv("CONSTANT_SEED", "11")

        first = self.framework().fit_atomic_system(
            self.water(), number_of_additions=5
        )
        second = self.framework().fit_atomic_system(
            self.water(), number_of_additions=5
        )

        for a, b in zip(first, second):
            assert np.array_equal(a.pos, b.pos)

    def test_iter_yields_each_addition(self, monkeypatch):
        monkeypatch.setenv("CONSTANT_SEED", "3")
        systems = self.framework().iter_fit_atomic_system(
            self.water(), number_of_additions=3
        )

        assert isinstance(next(systems), AtomicSystem)
        assert len(list(systems)) == 2
