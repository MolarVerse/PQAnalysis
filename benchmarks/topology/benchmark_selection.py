import numpy as np
import pytest

from PQAnalysis.core import Atom
from PQAnalysis.topology import Selection, Topology

ATOM_NAMES = ["O", "H", "H", "C", "N"]


def _topology(n_atoms, seed=42):
    rng = np.random.default_rng(seed)
    names = rng.choice(ATOM_NAMES, n_atoms)

    return Topology(atoms=[Atom(name) for name in names])


@pytest.mark.benchmark(group="Selection")
class BenchmarkSelection:

    @pytest.mark.parametrize(
        "selection", ["all", "O", "C|N", "H", "0..999"],
        ids=["all", "element", "union", "hydrogens", "range"],
    )
    @pytest.mark.parametrize("n_atoms", [3000, 30000])
    def benchmark_select(self, benchmark, selection, n_atoms):
        topology = _topology(n_atoms)

        benchmark(Selection(selection).select, topology)

    def benchmark_select_everything(self, benchmark):
        topology = _topology(30000)

        benchmark(Selection(None).select, topology)
