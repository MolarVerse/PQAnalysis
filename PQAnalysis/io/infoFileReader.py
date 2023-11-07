"""
A module containing the InfoFileReader class.

...

Classes
-------
InfoFileReader
    A class to read info files.
"""

from beartype.typing import Tuple, Dict

from .base import BaseReader


class InfoFileReader(BaseReader):
    """
    A class to read info files.

    Parameters
    ----------
    BaseReader : BaseReader
        A base class for all readers.

    Attributes
    ----------
    filename : str
        The name of the file to read from.
    format : str
        The format of the info file.
    """

    formats = ["pimd-qmcf", "qmcfc"]

    def __init__(self, filename: str, format: str = "pimd-qmcf") -> None:
        """
        Initializes the InfoFileReader with the given filename.

        Parameters
        ----------
        filename : str
            The name of the file to read from.
        format : str, optional
            The format of the info file, by default "pimd-qmcf"

        Raises
        ------
        ValueError
            If the format is not supported.
        """
        super().__init__(filename)

        if format not in self.formats:
            raise ValueError(
                f"Format {format} is not supported. Supported formats are {self.formats}.")

        self.format = format

    def read(self) -> Tuple[Dict, Dict | None]:
        """
        Reads the info file.

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
        if self.format == "pimd-qmcf":
            return self.read_pimd_qmcf()
        elif self.format == "qmcfc":
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
        """
        info = {}
        units = {}

        with open(self.filename, "r") as file:

            entry_counter = 0

            for line in file:
                line = line.split()

                if len(line) == 8:
                    info[line[1]] = entry_counter
                    units[line[1]] = line[3]
                    entry_counter += 1
                    info[line[4]] = entry_counter
                    units[line[4]] = line[6]
                    entry_counter += 1

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
        """
        info = {}

        with open(self.filename, "r") as file:

            entry_counter = 0

            for line in file:
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

        return info, None
