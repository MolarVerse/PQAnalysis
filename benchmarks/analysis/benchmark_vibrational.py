import numpy as np
import pytest

from PQAnalysis.analysis.vibrational.vibrational_analysis import calculate


def _molecule(n_atoms, seed=42):
    """A seeded positive semi-definite Hessian of a random molecule."""

    rng = np.random.default_rng(seed)
    coords = rng.uniform(-5.0, 5.0, (n_atoms, 3))
    masses = rng.uniform(1.0, 16.0, n_atoms)
    charges = rng.uniform(-0.5, 0.5, n_atoms)
    factor = rng.normal(size=(3 * n_atoms, 3 * n_atoms))

    return masses, coords, factor @ factor.T / (3 * n_atoms), charges


@pytest.mark.benchmark(group="Vibrational")
class BenchmarkVibrational:

    @pytest.mark.parametrize("n_atoms", [10, 30, 90])
    def benchmark_calculate(self, benchmark, n_atoms):
        masses, coords, hessian, _ = _molecule(n_atoms)

        benchmark(calculate, masses, coords, hessian, hessian_sign="positive")

    def benchmark_calculate_with_intensities(self, benchmark):
        masses, coords, hessian, charges = _molecule(30)

        benchmark(
            calculate,
            masses,
            coords,
            hessian,
            atom_charges=charges,
            hessian_sign="positive",
        )
