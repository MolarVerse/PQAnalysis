import numpy as np
import pytest

from PQAnalysis.analysis.spectrum_broadening import broaden, wavenumber_grid



@pytest.mark.benchmark(group="SpectrumBroadening")
class BenchmarkSpectrumBroadening:

    def benchmark_broaden_gaussian(self, benchmark):
        rng = np.random.default_rng(42)
        wavenumbers = rng.uniform(10.0, 4000.0, 200)
        intensities = rng.uniform(0.0, 10.0, 200)
        grid = wavenumber_grid()

        benchmark(broaden, wavenumbers, intensities, grid)

    def benchmark_broaden_lorentzian(self, benchmark):
        rng = np.random.default_rng(42)
        wavenumbers = rng.uniform(10.0, 4000.0, 200)
        intensities = rng.uniform(0.0, 10.0, 200)
        grid = wavenumber_grid()

        benchmark(
            broaden, wavenumbers, intensities, grid, kernel="lorentzian"
        )
