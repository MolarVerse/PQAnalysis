from ..core.atom import Atom


class ShakeTopology:
    def __init__(self, atoms: list[Atom] | None = None) -> None:
        if atoms is None:
            atoms = []

        self.atoms = atoms

    @property
    def atoms(self) -> list[Atom]:
        """
        The atoms in the system.

        Returns
        -------
        list[Atom]
            The atoms in the system.
        """
        return self._atoms

    @atoms.setter
    def atoms(self, atoms: list[Atom]) -> None:
        """
        The atoms in the system.

        Parameters
        ----------
        atoms : list[Atom]
            The atoms in the system.
        """
        self._atoms = atoms
