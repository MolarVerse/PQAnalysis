"""
A module containing the BondedTopology class. 

The BondedTopology class represents a bonded topology in a molecular
topology and is used to store the bonds, angles, dihedrals,
impropers, and shake bonds of a molecular system along with their properties.
"""

import logging

from beartype.typing import List

from PQAnalysis.types import PositiveInt
from PQAnalysis.exceptions import PQValueError
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

from .bond import Bond
from .angle import Angle
from .dihedral import Dihedral
from ._topology_properties import TopologyPropertiesMixin



class BondedTopology(TopologyPropertiesMixin):

    """
    A class to represent a bonded topology in a molecular topology.

    It inherits from TopologyPropertiesMixin.
    A mixin class to add the most common properties of a topology.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        bonds: List[Bond] | None = None,
        angles: List[Angle] | None = None,
        dihedrals: List[Dihedral] | None = None,
        impropers: List[Dihedral] | None = None,
        shake_bonds: List[Bond] | None = None,
        ordering_keys: List[str] | None = None
    ) -> None:
        """
        Parameters
        ----------
        bonds : List[Bond], optional
            A list of bonds in the bonded topology, by default None.
        angles : List[Angle], optional
            A list of angles in the bonded topology, by default None.
        dihedrals : List[Dihedral], optional
            A list of dihedrals in the bonded topology, by default None.
        impropers : List[Dihedral], optional
            A list of impropers in the bonded topology, by default None.
        shake_bonds : List[Bond], optional
            A list of shake bonds in the bonded topology, by default None.
        ordering_keys : List[str], optional
            A list of keys to order the blocks of the bonded topology,
            when initializing from a file, by default None.
        """

        self.bonds = bonds or []
        self.angles = angles or []
        self.dihedrals = dihedrals or []
        self.impropers = impropers or []
        self.shake_bonds = shake_bonds or []

        self.ordering_keys = ordering_keys

    @runtime_type_checking
    def extend_shake_bonds(
        self,
        shake_bonds: List[Bond],
        n_atoms: PositiveInt,
        n_extensions: PositiveInt = 1,
        n_atoms_per_extension: PositiveInt | None = None
    ) -> None:
        """
        Extend the shake bonds in the bonded topology.

        This function is useful when two or more atomic systems are merged,
        and the shake bonds need to be extended to the new atoms.

        Parameters
        ----------
        shake_bonds : List[Bond]
            A list of shake bonds to extend.
        n_atoms : PositiveInt
            The number of atoms in the system before the extension.
        n_extensions : PositiveInt, optional
            The number of extensions to perform, by default 1.
        n_atoms_per_extension : PositiveInt | None, optional
            The number of atoms in the system after each extension, by default None.

        Raises
        ------
        PQValueError
            If n_atoms_per_extension is not provided and n_extensions is not 1.
        PQValueError
            If n_atoms_per_extension is less than the highest index in the provided shake bonds.
        """

        if n_extensions != 1 and n_atoms_per_extension is None:
            self.logger.error(
                "n_atoms_per_extension must be provided if n_extensions is not 1.",
                exception=PQValueError
            )

        max_index_shake_bonds = max(
            [bond.index1
             for bond in shake_bonds] + [bond.index2 for bond in shake_bonds]
        )

        if n_atoms_per_extension is None:
            n_atoms_per_extension = 0
        elif n_atoms_per_extension < max_index_shake_bonds:
            self.logger.error(
                (
                    "n_atoms_per_extension must be greater or equal "
                    "than the highest index in the provided shake bonds."
                ),
                exception=PQValueError
            )

        for i in range(n_extensions):

            _n_atoms = n_atoms_per_extension * i + n_atoms

            for bond in shake_bonds:
                _bond = bond.copy()
                _bond.index1 += _n_atoms
                _bond.index2 += _n_atoms
                self.shake_bonds.append(_bond)
