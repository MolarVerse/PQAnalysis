"""
A module for writing topology files containing bonded information for the PQ and
QMCFC MD software packages. For more information about the topology file structure
please visit the documentation page of PQ https://molarverse.github.io/PQ/.
"""

import logging

from _io import TextIOWrapper as File  # type: ignore

from beartype.typing import List

from PQAnalysis.io.base import BaseWriter
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.topology import BondedTopology, Topology
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

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
        filename: str,
        mode: str | FileWritingMode = "w"
    ) -> None:
        """
        Parameters
        ----------
        filename : str
            The name of the topology file
        mode: str | FileWritingMode
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
            "shake": self._write_shake_info
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

        if isinstance(bonded_topology,
            Topology) and bonded_topology.bonded_topology is not None:
            bonded_topology = bonded_topology.bonded_topology
        elif not isinstance(bonded_topology, BondedTopology):
            self.logger.error(
                "Invalid bonded topology.",
                exception=TopologyFileError
            )

        self.open()

        if bonded_topology.ordering_keys is not None:
            keys = bonded_topology.ordering_keys
        else:
            keys = self.key_topology_map.keys()

        for key in keys:
            self.key_topology_map[key](bonded_topology, self.file)

        self.close()

    @classmethod
    def _write_bond_info(
        cls,
        bonded_topology: BondedTopology,
        file: File
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
        """
        if len(bonded_topology.bonds) != 0:
            lines = cls._get_bond_lines(bonded_topology)
            for line in lines:
                print(line, file=file)

    @classmethod
    def _write_angle_info(
        cls,
        bonded_topology: BondedTopology,
        file: File
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
        """
        if len(bonded_topology.angles) != 0:
            lines = cls._get_angle_lines(bonded_topology)
            for line in lines:
                print(line, file=file)

    @classmethod
    def _write_dihedral_info(
        cls,
        bonded_topology: BondedTopology,
        file: File
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
        """
        if len(bonded_topology.dihedrals) != 0:
            lines = cls._get_dihedral_lines(bonded_topology)
            for line in lines:
                print(line, file=file)

    @classmethod
    def _write_improper_info(
        cls,
        bonded_topology: BondedTopology,
        file: File
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
        cls,
        bonded_topology: BondedTopology,
        file: File
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
        n_linkers = len(bonded_topology.linkers)

        lines = []

        lines.append(
            f"BONDS {n_unique_indices} {n_unique_target_indices} {n_linkers}"
        )

        for bond in bonded_topology.bonds:
            lines.append(f"{bond.index1} {bond.index2} {bond.bond_type}")

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
        n_unique_indices3 = len(bonded_topology.unique_angle_target_indices)
        n_linkers = len(bonded_topology.angle_linkers)

        lines = []

        lines.append(
            f"ANGLES {n_unique_indices1} {n_unique_indices2} "
            f"{n_unique_indices3} {n_linkers}"
        )

        for angle in bonded_topology.angles:
            lines.append(
                f"{angle.index1} {angle.index2} "
                f"{angle.index3} {angle.angle_type}"
            )

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
            lines.append(
                f"{dihedral.index1} {dihedral.index2} {dihedral.index3} "
                f"{dihedral.index4} {dihedral.dihedral_type}"
            )

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
            lines.append(
                f"{improper.index1} {improper.index2} {improper.index3} "
                f"{improper.index4} {improper.improper_type}"
            )

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
            lines.append(
                f"{bond.index1} {bond.index2} "
                f"{bond.equilibrium_distance} {linker}"
            )

        lines.append("END")

        return lines
