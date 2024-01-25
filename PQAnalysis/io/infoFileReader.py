"""
A module containing the InfoFileReader class.
"""

from beartype.typing import Tuple, Dict

from .base import BaseReader
from PQAnalysis.traj import MDEngineFormat, MDEngineFormatError


class InfoFileReader(BaseReader):
    """
    This is a class to read info files from molecular dynamics simulations. The info file
    is a specific output file related to the energy file of PIMD-QMCF and QMCFC simulations.

    For more information on how the info file of PIMD-QMCF simulations is structured, see
    the corresponding documentation of the `PIMD-QMCF <https://molarverse.github.io/pimd_qmcf>`_ code.

    Calling the read method returns a tuple of dictionaries. For both dictionaries the keys are the names of the information strings (i.e. physical properties). The values of the first dictionary are the corresponding indices of the data entries, which can be used to index an :py:class:`~PQAnalysis.physicalData.energy.Energy` object. The values of the second dictionary are the corresponding units of the information strings (None if no units are given).
    """

    def __init__(self, filename: str, format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF) -> None:
        """
        Parameters
        ----------
        filename : str
            The name of the file to read from.
        format : MDEngineFormat | str, optional
            The format of the info file. Default is MDEngineFormat.PIMD_QMCF.
        """
        super().__init__(filename)

        self.format = MDEngineFormat(format)

    def read(self) -> Tuple[Dict, Dict | None]:
        """
        Calling the read method returns a tuple of dictionaries. For both dictionaries the keys are the names of the information strings (i.e. physical properties). The values of the first dictionary are the corresponding indices of the data entries, which can be used to index an :py:class:`~PQAnalysis.physicalData.energy.Energy` object. The values of the second dictionary are the corresponding units of the information strings (None if no units are given).

        Returns
        -------
        dict
            The information strings of the info file as a dictionary.
            The keys are the names of the information strings. The values are the
            corresponding data entry (columns in energy file).
        dict
            The units of the info file as a dictionary. The keys are the names of the
            information strings. The values are the corresponding units.
        """
        if self.format == MDEngineFormat.PIMD_QMCF:
            return self.read_pimd_qmcf()
        elif self.format == MDEngineFormat.QMCFC:
            return self.read_qmcfc()

    def read_pimd_qmcf(self) -> Tuple[Dict, Dict]:
        """
        Reads the info file in pimd-qmcf format.

        Returns
        -------
        dict
            The information strings of the info file as a dictionary.
            The keys are the names of the information strings. The values are the
            corresponding data entry (columns in energy file).
        dict
            The units of the info file as a dictionary. The keys are the names of the
            information strings. The values are the corresponding units.

        Raises
        ------
        MDEngineFormatError
            If the info file is not in pimd-qmcf format.
        """
        info = {}
        units = {}

        with open(self.filename, "r") as file:

            lines = file.readlines()

            entry_counter = 0

            for line in lines[3:-2]:
                line = line.split()

                if len(line) == 8:
                    info[line[1]] = entry_counter
                    units[line[1]] = line[3]
                    entry_counter += 1
                    info[line[4]] = entry_counter
                    units[line[4]] = line[6]
                    entry_counter += 1
                else:
                    raise MDEngineFormatError(
                        f"Info file {self.filename} is not in pimd-qmcf format.")

        return info, units

    def read_qmcfc(self) -> Tuple[Dict, None]:
        """
        Reads the info file in qmcfc format.

        Returns
        -------
        dict
            The information strings of the info file as a dictionary.
            The keys are the names of the information strings. The values are the
            corresponding data entry (columns in energy file).
        None
            For the qmcfc format, no units are given.

        Raises
        ------
        MDEngineFormatError
            If the info file is not in qmcfc format.
        """
        info = {}

        with open(self.filename, "r") as file:

            lines = file.readlines()

            entry_counter = 0

            for line in lines[3:-2]:
                line = line.split()

                if len(line) == 6:
                    info[line[1]] = entry_counter
                    entry_counter += 1
                    info[line[3]] = entry_counter
                    entry_counter += 1
                elif len(line) == 7:
                    info[' '.join(line[1:3])] = entry_counter
                    entry_counter += 1
                    info[line[4]] = entry_counter
                    entry_counter += 1
                else:
                    raise MDEngineFormatError(
                        f"Info file {self.filename} is not in qmcfc format.")

        return info, None
