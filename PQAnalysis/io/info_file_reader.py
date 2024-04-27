"""
A module containing the InfoFileReader class.
"""

from beartype.typing import Tuple, Dict

from PQAnalysis.traj import MDEngineFormat, MDEngineFormatError
from .base import BaseReader


class InfoFileReader(BaseReader):
    """
    This is a class to read info files from molecular dynamics simulations. The info file
    is a specific output file related to the energy file of PQ and QMCFC simulations.

    For more information on how the info file of PQ simulations is structured, see
    the corresponding documentation of the `PQ <https://molarverse.github.io/PQ>`_ code.

    Calling the read method returns a tuple of dictionaries. For both dictionaries
    the keys are the names of the information strings (i.e. physical properties). 
    The values of the first dictionary are the corresponding indices of the data
    entries, which can be used to index an :py:class:`~PQAnalysis.physicalData.energy.Energy`
    object. The values of the second dictionary are the corresponding units of the 
    information strings (None if no units are given).
    """

    def __init__(self,
                 filename: str,
                 engine_format: MDEngineFormat | str = MDEngineFormat.PQ
                 ) -> None:
        """
        Parameters
        ----------
        filename : str
            The name of the file to read from.
        engine_format : MDEngineFormat | str, optional
            The format of the info file. Default is MDEngineFormat.PQ.
        """
        super().__init__(filename)

        self.format = MDEngineFormat(engine_format)

    def read(self) -> Tuple[Dict, Dict | None]:
        """
        Calling the read method returns a tuple of dictionaries. For both dictionaries 
        the keys are the names of the information strings (i.e. physical properties).
        The values of the first dictionary are the corresponding indices of the data 
        entries, which can be used to index an :py:class:`~PQAnalysis.physicalData.energy.Energy`
        object. The values of the second dictionary are the corresponding units of the 
        information strings (None if no units are given).

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
            If the info file is not in PQ or qmcfc format.
        """
        if self.format == MDEngineFormat.PQ:
            return self.read_pq()

        if self.format == MDEngineFormat.QMCFC:
            return self.read_qmcfc()

        # should never reach this point - if it does, it is a bug
        raise MDEngineFormatError(
            f"Info file {self.filename} is not in PQ or qmcfc format."
        )

    def read_pq(self) -> Tuple[Dict, Dict]:
        """
        Reads the info file in PQ format.

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
            If the info file is not in PQ format.
        """
        info = {}
        units = {}

        with open(self.filename, "r", encoding='utf-8') as file:

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
                        f"Info file {self.filename} is not in PQ format.")

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

        with open(self.filename, "r", encoding='utf-8') as file:

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
