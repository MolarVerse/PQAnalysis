import pytest
import numpy as np

from PQAnalysis.core.cell import Cell



@pytest.mark.benchmark(group="Cell")
class BenchmarkCell:

    def benchmark_image_orthorhombic(self, benchmark):
        cell = Cell(10, 10, 10, 90, 90, 90)
        pos = np.array([5, 5, 5])
        benchmark(cell.image, pos)

    def benchmark_image_triclinic(self, benchmark):
        cell = Cell(10, 10, 10, 90.00000000001, 90.00000000001, 90.00000000001)
        pos = np.array([5, 5, 5])
        benchmark(cell.image, pos)
