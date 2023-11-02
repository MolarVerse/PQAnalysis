"""
A module containing the InfoFileReader class.

...

Classes
-------
InfoFileReader
    A class to read info files.
"""

from PQAnalysis.io.base import BaseReader


class InfoFileReader(BaseReader):
    """
    A class to read info files.

    Parameters
    ----------
    BaseReader : BaseReader
        A base class for all readers.
    """

    def __init__(self, filename: str):
        """
        Initializes the InfoFileReader with the given filename.

        Parameters
        ----------
        filename : str
            The name of the file to read from.
        """
        super().__init__(filename)

    def read(self) -> (dict, dict):
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
