from ..core.topology import Topology
from ..core.atomicSystem import AtomicSystem
from ..core.atom import Atom


class Frame:
    def __init__(self, system: AtomicSystem = AtomicSystem(), topology: Topology = None):
        self.system = system
        self.topology = topology

    @property
    def PBC(self):
        return self.system.PBC

    @property
    def cell(self):
        return self.system.cell

    @property
    def n_atoms(self):
        return self.system.n_atoms

    @property
    def pos(self):
        return self.system.pos

    def compute_com_frame(self, group=None):
        if group is None:
            group = self.n_atoms

        elif self.n_atoms % group != 0:
            raise ValueError(
                'Number of atoms in selection is not a multiple of group.')

        pos = np.zeros((0, 3))
        names = []

        j = 0
        for i in range(0, self.n_atoms, group):
            atomic_system = AtomicSystem(
                atoms=self.atoms[i:i+group], pos=self.pos[i:i+group], cell=self.cell)

            np.append(pos, atomic_system.compute_com(), axis=0)
            names.append(atomic_system.combined_name)

            j += 1

        names = [Atom(name, use_guess_type=False) for name in names]

        return Frame(AtomicSystem(pos, names, self.cell))
