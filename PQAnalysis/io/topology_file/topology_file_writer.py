"""
A module for writing topology files containing bonded information for the PQ and
QMCFC MD software packages. For more information about the topology file structure
please visit the documentation page of PQ https://molarverse.github.io/PQ/.
"""

import logging

from _io import TextIOWrapper as File  # type: ignore

from beartype.typing import Iterable, List

from PQAnalysis import __package_name__
from PQAnalysis.io.base import BaseWriter
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.topology import BondedTopology, Topology
from PQAnalysis.topology.bonded_topology import (
    Bond,
    Angle,
    Dihedral,
    JCoupling,
)

from .exceptions import TopologyFileError



class TopologyFileWriter(BaseWriter):

    """
    Class for writing topology files containing bonded information for the PQ and 
    QMCFC MD software packages. For more information about the topology file 
    structure please visit the documentation page of PQ https://molarverse.github.io/PQ/.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        filename: str | None = None,
        mode: str | FileWritingMode = "w",
    ) -> None:
        """
        Parameters
        ----------
        filename : str, optional
            The name of the topology file, Default is None, which means that the
            output is printed to stdout.
        mode: str | FileWritingMode, optional
            The writing mode. Default is "w".
            The writing mode can be either a string or a FileWritingMode enum value.
            Possible values are:
            - "w" or FileWritingMode.WRITE: write mode (default, no overwrite)
            - "a" or FileWritingMode.APPEND: append mode
            - "o" or FileWritingMode.OVERWRITE: overwrite mode
        """

        super().__init__(filename, mode=mode)

        self.key_topology_map = {
            "bonds": self._write_bond_info,
            "angles": self._write_angle_info,
            "dihedrals": self._write_dihedral_info,
            "impropers": self._write_improper_info,
            "shake": self._write_shake_info,
            "j_couplings": self._write_j_coupling_info,
            "dist_constraints": self._write_distance_constraint_info
        }

    @runtime_type_checking
    def write(self, bonded_topology: Topology | BondedTopology) -> None:
        """
        Writes the bonded topology to the file.

        Parameters
        ----------
        bonded_topology : Topology | BondedTopology
            The bonded topology to write to the file. If a Topology object is
            provided, the bonded topology will be extracted from it.

        Raises
        ------
        ValueError
            If the bonded topology is not a Topology or BondedTopology object.
        """

        if isinstance(
            bonded_topology, Topology
        ) and bonded_topology.bonded_topology is not None:
            bonded_topology = bonded_topology.bonded_topology
        elif not isinstance(bonded_topology, BondedTopology):
            self.logger.error(
                "Invalid bonded topology.", exception=TopologyFileError
            )

        if bonded_topology.ordering_keys is not None:
            keys = bonded_topology.ordering_keys
        else:
            keys = self.key_topology_map.keys()

        self._check_types_given(bonded_topology, keys)

        self.open()

        try:
            for key in keys:
                self.key_topology_map[key](bonded_topology, self.file)
        finally:
            self.close()

    @classmethod
    def _check_types_given(
        cls, bonded_topology: BondedTopology, keys: Iterable[str]
    ) -> None:
        """
        Checks that a type is defined for all bonded parameters to be written.

        This check has to be performed before the output file is opened,
        because opening it truncates an already existing file.

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object to check.
        keys : Iterable[str]
            The keys of the topology sections that will be written.

        Raises
        ------
        TopologyFileError
            If any bond, angle, dihedral, improper or j-coupling to be
            written does not have a type defined.
        """

        type_names = {
            "bonds": "bond",
            "angles": "angle",
            "dihedrals": "dihedral",
            "impropers": "improper",
            "j_couplings": "j-coupling",
        }

        for key in keys:
            if key in type_names:
                cls._check_type_given(
                    getattr(bonded_topology, key), type_names[key]
                )

    @classmethod
    def _write_bond_info(
        cls, bonded_topology: BondedTopology, file: File
    ) -> None:
        """
        Determines if the bonded topology contains bonds and 
        writes the bond information to the file.

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the bond information.
        file : File
            The file object to write the bond information to.
            
        Raises
        ------
        TopologyFileError
            If any bond in the bonded topology does not have a bond type defined.
        """

        if len(bonded_topology.bonds) != 0:
            lines = cls._get_bond_lines(bonded_topology)
            for line in lines:
                print(line, file=file)

    @classmethod
    def _write_angle_info(
        cls, bonded_topology: BondedTopology, file: File
    ) -> None:
        """
        Determines if the bonded topology contains angles and
        writes the angle information to the file.

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the angle information.
        file : File 
            The file object to write the angle information to.
        
        Raises
        ------
        TopologyFileError
            If any angle in the bonded topology does not have an angle type defined.
        """

        if len(bonded_topology.angles) != 0:
            lines = cls._get_angle_lines(bonded_topology)
            for line in lines:
                print(line, file=file)

    @classmethod
    def _write_dihedral_info(
        cls, bonded_topology: BondedTopology, file: File
    ) -> None:
        """
        Determines if the bonded topology contains dihedrals and
        writes the dihedral information to the file.

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the dihedral information.
        file : File
            The file object to write the dihedral information to.
            
        Raises
        ------
        TopologyFileError
            If any dihedral in the bonded topology does not have a dihedral type defined.
        """

        if len(bonded_topology.dihedrals) != 0:
            lines = cls._get_dihedral_lines(bonded_topology)
            for line in lines:
                print(line, file=file)

    @classmethod
    def _write_improper_info(
        cls, bonded_topology: BondedTopology, file: File
    ) -> None:
        """
        Writes the improper information to the file.

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the improper information.
        file : File
            The file object to write the improper information to.
        """

        if len(bonded_topology.impropers) != 0:
            lines = cls._get_improper_lines(bonded_topology)
            for line in lines:
                print(line, file=file)

    @classmethod
    def _write_shake_info(
        cls, bonded_topology: BondedTopology, file: File
    ) -> None:
        """
        Writes the shake information to the file.

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the shake information.
        file : File
            The file object to write the shake information to.
        """
        if len(bonded_topology.shake_bonds) != 0:
            lines = cls._get_shake_lines(bonded_topology)
            for line in lines:
                print(line, file=file)

    @classmethod
    def _write_j_coupling_info(
        cls, bonded_topology: BondedTopology, file: File
    ) -> None:
        """
        Determines if the bonded topology contains j-couplings and
        writes the j-coupling information to the file.

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the j-coupling information.
        file : File
            The file object to write the j-coupling information to.

        Raises
        ------
        TopologyFileError
            If any j-coupling in the bonded topology does not have
            a j-coupling type defined.
        """

        if len(bonded_topology.j_couplings) != 0:
            lines = cls._get_j_coupling_lines(bonded_topology)
            for line in lines:
                print(line, file=file)

    @classmethod
    def _write_distance_constraint_info(
        cls, bonded_topology: BondedTopology, file: File
    ) -> None:
        """
        Writes the distance constraint information to the file.

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the distance
            constraint information.
        file : File
            The file object to write the distance constraint information to.
        """
        if len(bonded_topology.distance_constraints) != 0:
            lines = cls._get_distance_constraint_lines(bonded_topology)
            for line in lines:
                print(line, file=file)

    @classmethod
    def _check_type_given(
        cls, types: List[Bond | Angle | Dihedral | JCoupling], type_name: str
    ) -> None:
        """
        Check if the type is given for each bond, angle, dihedral, or j-coupling.

        Parameters
        ----------
        types : List[Bond | Angle | Dihedral | JCoupling]
            The list of bonds, angles, dihedrals, or j-couplings to check.
        type_name : str
            The name of the type to check.
        """

        def get_type(type_: Bond | Angle | Dihedral | JCoupling) -> int | None:

            if isinstance(type_, Bond):
                return type_.bond_type

            if isinstance(type_, Angle):
                return type_.angle_type

            if isinstance(type_, JCoupling):
                return type_.j_coupling_type

            return type_.dihedral_type

        if any(get_type(type_) is None for type_ in types):
            cls.logger.error(
                (
                    f"In order to write the {type_name} information in 'PQ' topology format, "
                    f"all {type_name}s must have a {type_name} type defined."
                ),
                exception=TopologyFileError
            )

    @staticmethod
    def _get_bond_lines(bonded_topology: BondedTopology) -> List[str]:
        """
        Get the bond lines for the bonded topology.

        The lines contain one header line, one line for each bond,
        and an end line in the following format:

        BONDS n_unique_indices n_unique_target_indices n_linkers
        index1 index2 bond_type
        ...
        END

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the bond information.

        Returns
        -------
        List[str]
            The list of bond lines.
        """
        n_unique_indices = len(bonded_topology.unique_bond1_indices)
        n_unique_target_indices = len(bonded_topology.unique_bond2_indices)
        n_linkers = len(bonded_topology.bond_linkers)

        lines = []

        lines.append(
            f"BONDS {n_unique_indices} {n_unique_target_indices} {n_linkers}"
        )

        for bond in bonded_topology.bonds:
            line = f"{bond.index1:>5d} {bond.index2:>5d} {bond.bond_type:>5d}"

            if bond.is_linker:
                line += " *"

            if bond.comment is not None:
                line += f" # {bond.comment}"

            lines.append(line)

        lines.append("END")

        return lines

    @staticmethod
    def _get_angle_lines(bonded_topology: BondedTopology) -> List[str]:
        """
        Get the angle lines for the bonded topology.

        The lines contain one header line, one line for each angle,
        and an end line in the following format:

        ANGLES n_unique_indices1 n_unique_indices2 n_unique_indices3 n_linkers
        index1 index2 index3 angle_type
        ...
        END

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the angle information.

        Returns
        -------
        List[str]
            The list of angle lines.
        """
        n_unique_indices1 = len(bonded_topology.unique_angle1_indices)
        n_unique_indices2 = len(bonded_topology.unique_angle2_indices)
        n_unique_indices3 = len(bonded_topology.unique_angle3_indices)
        n_linkers = len(bonded_topology.angle_linkers)

        lines = []

        lines.append(
            f"ANGLES {n_unique_indices1} {n_unique_indices2} "
            f"{n_unique_indices3} {n_linkers}"
        )

        for angle in bonded_topology.angles:
            line = (
                f"{angle.index1:>5d} {angle.index2:>5d} "
                f"{angle.index3:>5d} {angle.angle_type:>5d}"
            )

            if angle.is_linker:
                line += " *"

            if angle.comment is not None:
                line += f" # {angle.comment}"

            lines.append(line)

        lines.append("END")

        return lines

    @staticmethod
    def _get_dihedral_lines(bonded_topology: BondedTopology) -> List[str]:
        """
        Get the dihedral lines for the bonded topology.

        The lines contain one header line, one line for each dihedral,
        and an end line in the following format:

        DIHEDRALS n_unique_indices1 n_unique_indices2 n_unique_indices3 n_unique_indices4
        index1 index2 index3 index4 dihedral_type
        ...
        END

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the dihedral information.

        Returns
        -------
        List[str]
            The list of dihedral lines.
        """
        n_unique_indices1 = len(bonded_topology.unique_dihedral1_indices)
        n_unique_indices2 = len(bonded_topology.unique_dihedral2_indices)
        n_unique_indices3 = len(bonded_topology.unique_dihedral3_indices)
        n_unique_indices4 = len(bonded_topology.unique_dihedral4_indices)

        lines = []

        lines.append(
            f"DIHEDRALS {n_unique_indices1} {n_unique_indices2} "
            f"{n_unique_indices3} {n_unique_indices4}"
        )

        for dihedral in bonded_topology.dihedrals:
            line = (
                f"{dihedral.index1:>5d} {dihedral.index2:>5d} "
                f"{dihedral.index3:>5d} {dihedral.index4:>5d} "
                f"{dihedral.dihedral_type:>5d}"
            )

            if dihedral.is_linker:
                line += " *"

            if dihedral.comment is not None:
                line += f" # {dihedral.comment}"

            lines.append(line)

        lines.append("END")

        return lines

    @staticmethod
    def _get_improper_lines(bonded_topology: BondedTopology) -> List[str]:
        """
        Get the improper lines for the bonded topology.

        The lines contain one header line, one line for each improper,
        and an end line in the following format:

        IMPROPERS n_unique_indices1 n_unique_indices2 n_unique_indices3 n_unique_indices4
        index1 index2 index3 index4 improper_type
        ...
        END

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the improper information.

        Returns
        -------
        List[str]
            The list of improper lines.
        """
        n_unique_indices1 = len(bonded_topology.unique_improper1_indices)
        n_unique_indices2 = len(bonded_topology.unique_improper2_indices)
        n_unique_indices3 = len(bonded_topology.unique_improper3_indices)
        n_unique_indices4 = len(bonded_topology.unique_improper4_indices)

        lines = []

        lines.append(
            f"IMPROPERS {n_unique_indices1} {n_unique_indices2} "
            f"{n_unique_indices3} {n_unique_indices4}"
        )

        for improper in bonded_topology.impropers:

            line = (
                f"{improper.index1:>5d} {improper.index2:>5d} "
                f"{improper.index3:>5d} {improper.index4:>5d} "
                f"{improper.dihedral_type:>5d}"
            )

            if improper.is_linker:
                line += " *"

            if improper.comment is not None:
                line += f" # {improper.comment}"

            lines.append(line)

        lines.append("END")

        return lines

    @staticmethod
    def _get_shake_lines(bonded_topology: BondedTopology) -> List[str]:
        """
        Get the shake lines for the bonded topology.

        The lines contain one header line, one line for each shake bond,
        and an end line in the following format:

        SHAKE n_unique_indices n_unique_target_indices n_linkers
        index1 index2 equilibrium_distance linker
        ...
        END

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the shake information.

        Returns
        -------
        List[str]
            The list of shake lines.
        """
        n_unique_indices = len(bonded_topology.unique_shake_indices)
        n_unique_target_indices = len(
            bonded_topology.unique_shake_target_indices
        )
        n_linkers = len(bonded_topology.shake_linkers)

        lines = []

        lines.append(
            f"SHAKE {n_unique_indices} {n_unique_target_indices} {n_linkers}"
        )

        for bond in bonded_topology.shake_bonds:
            linker = "*" if bond.is_linker else ""

            line = (
                f"{bond.index1:>5d} {bond.index2:>5d} "
                f"{bond.equilibrium_distance:16.12f}\t{linker}"
            )

            if bond.comment is not None:
                line += f" # {bond.comment}"

            lines.append(line)

        lines.append("END")

        return lines

    @staticmethod
    def _get_j_coupling_lines(bonded_topology: BondedTopology) -> List[str]:
        """
        Get the j-coupling lines for the bonded topology.

        The lines contain one header line, one line for each j-coupling,
        and an end line in the following format:

        J_COUPLINGS n_unique_indices1 n_unique_indices2 n_unique_indices3 n_unique_indices4
        index1 index2 index3 index4 j_coupling_type
        ...
        END

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the j-coupling information.

        Returns
        -------
        List[str]
            The list of j-coupling lines.
        """
        n_unique_indices1 = len(bonded_topology.unique_j_coupling1_indices)
        n_unique_indices2 = len(bonded_topology.unique_j_coupling2_indices)
        n_unique_indices3 = len(bonded_topology.unique_j_coupling3_indices)
        n_unique_indices4 = len(bonded_topology.unique_j_coupling4_indices)

        lines = []

        lines.append(
            f"J_COUPLINGS {n_unique_indices1} {n_unique_indices2} "
            f"{n_unique_indices3} {n_unique_indices4}"
        )

        for j_coupling in bonded_topology.j_couplings:
            line = (
                f"{j_coupling.index1:>5d} {j_coupling.index2:>5d} "
                f"{j_coupling.index3:>5d} {j_coupling.index4:>5d} "
                f"{j_coupling.j_coupling_type:>5d}"
            )

            if j_coupling.comment is not None:
                line += f" # {j_coupling.comment}"

            lines.append(line)

        lines.append("END")

        return lines

    @staticmethod
    def _get_distance_constraint_lines(
        bonded_topology: BondedTopology
    ) -> List[str]:
        """
        Get the distance constraint lines for the bonded topology.

        The lines contain one header line, one line for each distance
        constraint, and an end line in the following format:

        DIST_CONSTRAINTS n_unique_indices n_unique_target_indices
        index1 index2 lower_distance upper_distance spring_constant dk/dt
        ...
        END

        Parameters
        ----------
        bonded_topology : BondedTopology
            The bonded topology object containing the distance
            constraint information.

        Returns
        -------
        List[str]
            The list of distance constraint lines.
        """
        n_unique_indices = len(
            bonded_topology.unique_distance_constraint1_indices
        )
        n_unique_target_indices = len(
            bonded_topology.unique_distance_constraint2_indices
        )

        lines = []

        lines.append(
            f"DIST_CONSTRAINTS {n_unique_indices} {n_unique_target_indices}"
        )

        for constraint in bonded_topology.distance_constraints:
            line = (
                f"{constraint.index1:>5d} {constraint.index2:>5d} "
                f"{constraint.lower_distance:16.12f} "
                f"{constraint.upper_distance:16.12f} "
                f"{constraint.spring_constant:16.12f} "
                f"{constraint.d_spring_constant_dt:16.12f}"
            )

            if constraint.comment is not None:
                line += f" # {constraint.comment}"

            lines.append(line)

        lines.append("END")

        return lines
