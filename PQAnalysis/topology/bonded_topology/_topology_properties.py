"""
A module containing a mixin class to provide properties for a BondedTopology object.
"""

from beartype.typing import List, Set

from PQAnalysis.types import PositiveInt
from .bond import Bond
from .angle import Angle
from .dihedral import Dihedral



class TopologyPropertiesMixin:

    """
    A mixin class to add the most common properties of a topology.
    """

    @property
    def unique_bond1_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the first atoms in the bonds."""
        return {bond.index1 for bond in self.bonds}

    @property
    def unique_bond2_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the second atoms in the bonds."""
        return {bond.index2 for bond in self.bonds}

    @property
    def bond_linkers(self) -> List[Bond]:
        """List[Bond]: The bonds that are linkers."""
        return [bond for bond in self.bonds if bond.is_linker]

    @property
    def unique_angle1_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the first atoms in the angles."""
        return {angle.index1 for angle in self.angles}

    @property
    def unique_angle2_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the second atoms in the angles."""
        return {angle.index2 for angle in self.angles}

    @property
    def unique_angle3_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the third atoms in the angles."""
        return {angle.index3 for angle in self.angles}

    @property
    def angle_linkers(self) -> List[Angle]:
        """List[Angle]: The angles that are linkers."""
        return [angle for angle in self.angles if angle.is_linker]

    @property
    def unique_dihedral1_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the first atoms in the dihedrals."""
        return {dihedral.index1 for dihedral in self.dihedrals}

    @property
    def unique_dihedral2_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the second atoms in the dihedrals."""
        return {dihedral.index2 for dihedral in self.dihedrals}

    @property
    def unique_dihedral3_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the third atoms in the dihedrals."""
        return {dihedral.index3 for dihedral in self.dihedrals}

    @property
    def unique_dihedral4_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the fourth atoms in the dihedrals."""
        return {dihedral.index4 for dihedral in self.dihedrals}

    @property
    def dihedral_linkers(self) -> List[Dihedral]:
        """List[Dihedral]: The dihedrals that are linkers."""
        return [dihedral for dihedral in self.dihedrals if dihedral.is_linker]

    @property
    def unique_improper1_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the first atoms in the impropers."""
        return {improper.index1 for improper in self.impropers}

    @property
    def unique_improper2_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the second atoms in the impropers."""
        return {improper.index2 for improper in self.impropers}

    @property
    def unique_improper3_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the third atoms in the impropers."""
        return {improper.index3 for improper in self.impropers}

    @property
    def unique_improper4_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the fourth atoms in the impropers."""
        return {improper.index4 for improper in self.impropers}

    @property
    def improper_linkers(self) -> List[Dihedral]:
        """List[Dihedral]: The impropers that are linkers."""
        return [improper for improper in self.impropers if improper.is_linker]

    @property
    def unique_shake_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the atoms in the shake bonds."""
        return {bond.index1 for bond in self.shake_bonds}

    @property
    def unique_shake_target_indices(self) -> Set[PositiveInt]:
        """Set[PositiveInt]: The unique indices of the target atoms in the shake bonds."""
        return {bond.index2 for bond in self.shake_bonds}

    @property
    def shake_linkers(self) -> List[Bond]:
        """List[Bond]: The shake bonds that are linkers."""
        return [bond for bond in self.shake_bonds if bond.is_linker]
